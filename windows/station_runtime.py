"""Probe and stop only the local service that belongs to an installed station."""
from pathlib import Path
import argparse
import json
import socket
import time
import urllib.error
import urllib.request


SERVICE = '16s-workstation'


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, request, response, code, message, headers, new_url):
        return None


# Local service management must never pass through an HTTP proxy.
HTTP = urllib.request.build_opener(urllib.request.ProxyHandler({}), _NoRedirect())


def _port(value):
    if isinstance(value, bool) or not isinstance(value, int) or not 1024 <= value <= 65535:
        raise ValueError('端口必须为 1024–65535')
    return value


def _owned(root, health):
    if not isinstance(health, dict):
        return False
    try:
        program = Path(health['program_dir']).resolve()
        return (health.get('service_id') == SERVICE
                and Path(health['install_root']).resolve() == root
                and program.parent == root/'程序'
                and isinstance(health.get('pid'), int) and not isinstance(health['pid'], bool)
                and health['pid'] > 1
                and isinstance(health.get('instance_id'), str) and bool(health['instance_id'])
                and isinstance(health.get('busy'), bool))
    except (KeyError, TypeError, ValueError, OSError):
        return False


def probe(root, port):
    root, port = Path(root).resolve(), _port(port)
    try:
        request = urllib.request.Request(f'http://127.0.0.1:{port}/api/health')
        with HTTP.open(request, timeout=3) as response:
            health = json.loads(response.read(65537))
        if _owned(root, health):
            return dict(state='running', health=health)
        return dict(state='occupied', health={})
    except (urllib.error.URLError, TimeoutError, OSError, ValueError):
        with socket.socket() as connection:
            connection.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            try:
                connection.bind(('127.0.0.1', port))
            except OSError:
                return dict(state='occupied', health={})
        return dict(state='offline', health={})


def wait(root, port, seconds=90):
    end = time.monotonic()+seconds
    while True:
        result = probe(root, port)
        if result['state'] == 'running':
            return result
        if result['state'] == 'occupied':
            raise ValueError('端口已有其他程序；未打开它，请更换工作站端口。')
        if time.monotonic() >= end:
            raise ValueError('工作站启动超时，请查看启动终端和日志/server.log。')
        time.sleep(.3)


def stop(root, port):
    root, port = Path(root).resolve(), _port(port)
    value = probe(root, port)
    if value['state'] == 'offline':
        return value
    if value['state'] != 'running':
        raise ValueError('端口由其他程序使用，未停止它。请改用其他端口。')
    health = value['health']
    if health['busy']:
        raise ValueError('还有未结束的分析，请先完成或在网页中取消，未停止工作站。')
    record_file = root/'配置/service.json'
    try:
        if record_file.is_symlink():
            raise ValueError('服务记录不能是符号链接')
        record = json.loads(record_file.read_text(encoding='utf-8-sig'))
        if (not isinstance(record, dict) or not isinstance(record.get('token'), str)
                or not record['token'] or record.get('port') != port
                or record.get('pid') != health['pid']
                or record.get('instance_id') != health['instance_id']
                or Path(record['program_dir']).resolve() != Path(health['program_dir']).resolve()
                or Path(record['install_root']).resolve() != root):
            raise ValueError('服务记录不匹配')
    except (OSError, ValueError, KeyError, TypeError) as error:
        raise ValueError('服务进程与记录不一致，未停止它。请在原终端按 Ctrl+C。') from error
    request = urllib.request.Request(
        f'http://127.0.0.1:{port}/api/station-stop', data=b'{}', method='POST',
        headers={'X-Workstation-Token': record['token'], 'Content-Type': 'application/json'})
    try:
        with HTTP.open(request, timeout=5) as response:
            response.read(65537)
    except urllib.error.HTTPError as error:
        if error.code == 409:
            raise ValueError('还有未结束的分析，请先完成或在网页中取消，未停止工作站。') from error
        raise ValueError('服务拒绝停止请求；未强制终止，请查看原终端。') from error
    except (urllib.error.URLError, TimeoutError, OSError) as error:
        raise ValueError('无法确认停止请求；未强制终止，请查看原终端。') from error
    end = time.monotonic()+15
    while time.monotonic() < end:
        current = probe(root, port)
        if current['state'] == 'offline':
            return current
        if (current['state'] == 'running'
                and current['health']['instance_id'] != health['instance_id']):
            raise ValueError('端口服务已改变，未停止新服务，请查看原终端。')
        time.sleep(.1)
    raise ValueError('服务尚未退出，未强制终止，请查看原终端。')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['probe', 'wait', 'stop'])
    parser.add_argument('root')
    parser.add_argument('port', type=int)
    parser.add_argument('--seconds', type=int, default=90)
    args = parser.parse_args()
    try:
        result = (wait(args.root, args.port, args.seconds) if args.action == 'wait'
                  else stop(args.root, args.port) if args.action == 'stop'
                  else probe(args.root, args.port))
        print(json.dumps(result, ensure_ascii=False))
    except (ValueError, OSError) as error:
        parser.exit(1, str(error)+'\n')
