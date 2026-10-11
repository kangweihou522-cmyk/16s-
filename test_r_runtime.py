"""Check process boundaries around isolated R; fake Conda does not validate R packages."""
import importlib.util
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


SPEC = importlib.util.spec_from_file_location("workstation_r_runtime", Path(__file__).with_name("app.py"))
app = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(app)


class IsolatedRRuntimeTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="16s-r-isolation-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.conda = self.root / "fake conda"
        r_bin = self.root / "r bin"
        r_bin.mkdir()
        fake_r = r_bin / "Rscript"
        fake_r.write_text(
            f"#!{sys.executable}\n"
            "import json, os, sys\n"
            "keys = ['R_LIBS', 'R_LIBS_USER', 'R_LIBS_SITE', 'R_PROFILE_USER', 'R_ENVIRON_USER', 'R_HOME', 'LD_LIBRARY_PATH', 'LD_PRELOAD', 'PKG_CONFIG_PATH', 'TEST_KEEP_ME']\n"
            "print(json.dumps({'args': ['Rscript'] + sys.argv[1:], 'selected': os.environ.get('TEST_SELECTED_CONDA_ENV'), 'environment': {k: os.environ.get(k) for k in keys}}, ensure_ascii=False))\n",
            encoding="utf8",
        )
        fake_r.chmod(0o755)
        # A real process captures precisely the argv and environment it receives.
        # No R interpreter, packages, dependency installation, or network is faked
        # as a successful scientific analysis.
        self.conda.write_text(
            f"#!{sys.executable}\n"
            "import json, os, subprocess, sys\n"
            "args = sys.argv[1:]\n"
            "if args[:2] == ['env', 'list']:\n"
            "    print(json.dumps({'envs': ['/mock/envs/qiime-test', '/mock/envs/16s-workstation-r']}))\n"
            "    sys.exit(0)\n"
            "if '-n' not in args:\n"
            "    print('Named environment required', file=sys.stderr)\n"
            "    sys.exit(9)\n"
            "selected = args[args.index('-n') + 1]\n"
            "if selected == 'missing-r':\n"
            "    print('EnvironmentNameNotFound: missing-r', file=sys.stderr)\n"
            "    sys.exit(17)\n"
            "if 'Rscript' in args:\n"
            "    execution_env = os.environ.copy()\n"
            "    execution_env['TEST_SELECTED_CONDA_ENV'] = selected\n"
            # Conda activation can introduce defaults after Python has sanitized
            # its environment; validate the actual post-activation R process.
            "    execution_env.update({'R_LIBS': '/tmp/activation/library', 'R_LIBS_USER': '/tmp/activation/user-library', 'R_LIBS_SITE': '/tmp/activation/site-library', 'R_HOME': '/tmp/activation/wrong-R', 'LD_LIBRARY_PATH': '/tmp/activation/wrong-lib', 'R_PROFILE_USER': '/tmp/activation/profile.R', 'R_ENVIRON_USER': '/tmp/activation/Renviron'})\n"
            f"    execution_env['PATH'] = {str(r_bin)!r} + os.pathsep + execution_env.get('PATH', '')\n"
            "    sys.exit(subprocess.run(args[args.index('-n') + 2:], env=execution_env).returncode)\n"
            "else:\n"
            "    print(json.dumps({'args': args, 'selected': selected, 'environment': {'R_LIBS_USER': os.environ.get('R_LIBS_USER')}}))\n",
            encoding="utf8",
        )
        self.conda.chmod(0o755)
        self.config = {"conda": str(self.conda), "env": "qiime-test", "r_env": "isolated-r"}

    def test_r_has_own_default_environment_without_qiime_selection(self):
        config = {"conda": str(self.conda), "env": ""}
        args = app.r_command(config, "--version")
        self.assertEqual(args[args.index("-n") + 1], "16s-workstation-r")
        self.assertEqual(args[-2:], ["Rscript", "--version"])
        self.assertNotIn("--vanilla", args)

    def test_r_and_qiime_select_different_named_environments(self):
        r_args = app.r_command(self.config, "check_r.R", "DejaVu Serif")
        qiime_args = app.conda_cmd(self.config) + ["qiime", "info"]
        self.assertEqual(r_args[r_args.index("-n") + 1], "isolated-r")
        self.assertEqual(qiime_args[qiime_args.index("-n") + 1], "qiime-test")
        self.assertIn("--vanilla", r_args)
        self.assertNotIn("Rscript", qiime_args)
        self.assertEqual(self.config["env"], "qiime-test")

    def test_r_environment_name_rejects_paths_and_option_injection(self):
        for env in ("../qiime-test", "/tmp/r-env", "--prefix", "my env", "env\nname"):
            with self.subTest(env=env), self.assertRaises(ValueError):
                app.r_command(self.config | {"r_env": env}, "--version")

    def test_real_qiime_process_keeps_its_environment_and_receives_no_r_overrides(self):
        job = {"log": str(self.root / "qiime.log"), "run": str(self.root)}
        with patch.dict(os.environ, {"R_LIBS_USER": "/mock/qiime-own-r-library"}):
            app.q(job, self.config, "info")
        captured = next(json.loads(line) for line in Path(job["log"]).read_text(encoding="utf8").splitlines()
                        if line.startswith('{"args":'))
        self.assertEqual(captured["selected"], "qiime-test")
        self.assertEqual(captured["environment"]["R_LIBS_USER"], "/mock/qiime-own-r-library")
        self.assertEqual(captured["args"][-2:], ["qiime", "info"])
        self.assertNotIn("Rscript", captured["args"])
        self.assertNotIn("--vanilla", captured["args"])

    def test_real_subprocess_excludes_user_r_library_and_startup_overrides(self):
        polluted = {
            "R_LIBS": "/tmp/system-R-4.4/library",
            "R_LIBS_USER": "/tmp/user-R-4.4/library",
            "R_LIBS_SITE": "/tmp/site-R-4.4/library",
            "R_PROFILE_USER": "/tmp/untrusted-profile.R",
            "R_ENVIRON_USER": "/tmp/untrusted-Renviron",
            "R_HOME": "/tmp/system-R-4.4/home",
            "LD_LIBRARY_PATH": "/tmp/mixed-library",
            "LD_PRELOAD": "/tmp/nonexistent-preload.so",
            "PKG_CONFIG_PATH": "/tmp/mixed-pkgconfig",
            "TEST_KEEP_ME": "preserved unrelated variable",
        }
        script = self.root / "图 16S '脚本'.R"
        result_path = self.root / "result with spaces"
        result_path.mkdir()
        job = {"log": str(self.root / "workflow.log")}
        with patch.dict(os.environ, polluted):
            process_env = app.r_execution_env(self.config)
            app.command(job, app.r_command(self.config, script, result_path, "Font 'quoted'"),
                        cwd=result_path, process_env=process_env)
            # Preparing an isolated process must not contaminate the parent or
            # subsequent QIIME process environment.
            self.assertEqual(os.environ["R_LIBS_USER"], polluted["R_LIBS_USER"])
        captured = next(json.loads(line) for line in Path(job["log"]).read_text(encoding="utf8").splitlines()
                        if line.startswith('{"args":'))
        self.assertEqual(captured["selected"], "isolated-r")
        self.assertEqual(captured["args"][-3:], [str(script), str(result_path), "Font 'quoted'"])
        self.assertIn("--vanilla", captured["args"])
        for name in ("R_LIBS", "R_LIBS_USER", "R_LIBS_SITE"):
            self.assertEqual(captured["environment"][name], os.devnull, name)
        for name in ("R_PROFILE_USER", "R_ENVIRON_USER"):
            self.assertEqual(captured["environment"][name], os.devnull, name)
        self.assertEqual(captured["environment"]["TEST_KEEP_ME"], polluted["TEST_KEEP_ME"])
        for name in ("R_HOME", "LD_LIBRARY_PATH", "LD_PRELOAD", "PKG_CONFIG_PATH"):
            self.assertIsNone(captured["environment"][name], name)
        command_log = json.loads((self.root / "commands.jsonl").read_text(encoding="utf8"))
        self.assertEqual(command_log["exit_code"], 0)
        self.assertEqual(command_log["args"][-3:], [str(script), str(result_path), "Font 'quoted'"])
        self.assertNotIn("process_env", command_log)

    def test_missing_r_environment_fails_without_using_system_r(self):
        fallback_marker = self.root / "system-r-was-used"
        system_r = self.root / "Rscript"
        system_r.write_text(f"#!{sys.executable}\nfrom pathlib import Path\nPath({str(fallback_marker)!r}).touch()\n")
        system_r.chmod(0o755)
        config = self.config | {"r_env": "missing-r"}
        job = {"log": str(self.root / "missing.log")}
        with patch.dict(os.environ, {"PATH": str(self.root) + os.pathsep + os.environ.get("PATH", "")}):
            with self.assertRaisesRegex(RuntimeError, "17"):
                app.command(job, app.r_command(config, "--version"), process_env=app.r_execution_env(config))
        self.assertFalse(fallback_marker.exists())
        log = Path(job["log"]).read_text(encoding="utf8")
        self.assertIn("EnvironmentNameNotFound", log)
        command_log = json.loads((self.root / "commands.jsonl").read_text(encoding="utf8"))
        self.assertEqual(command_log["exit_code"], 17)
        self.assertEqual(command_log["args"][0], str(self.conda))

    def test_preflight_reports_missing_r_without_blocking_qiime_diagnostics(self):
        job = {"log": str(self.root / "preflight.log")}
        app.preflight(job, self.config | {"r_env": "missing-r", "font": "DejaVu Serif"})
        checks = {check["component"]: check for check in job["result"]["checks"]}
        self.assertTrue(checks["conda"]["ok"])
        self.assertTrue(checks["QIIME 2"]["ok"])
        self.assertFalse(checks["专用 R 环境与绘图依赖"]["ok"])
        self.assertIn("EnvironmentNameNotFound", checks["专用 R 环境与绘图依赖"]["details"])
        self.assertTrue(all(check["ok"] for name, check in checks.items() if name.startswith("plugin ")))


if __name__ == "__main__":
    unittest.main()
