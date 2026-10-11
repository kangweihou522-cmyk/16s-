"""Exercise the real launcher with isolated local Git remotes and tiny apps."""
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


@unittest.skipUnless(shutil.which("git"), "Git is required for updater integration tests")
class LauncherIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="16s-launch-test-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.remote = self.root / "remote.git"
        self.seed = self.root / "seed"
        self.checkout = self.root / "workstation"
        self.env = os.environ.copy()
        self.env.pop("WORKSTATION_UPDATED_ONCE", None)
        self.env.update({
            "GIT_CONFIG_NOSYSTEM": "1",
            "GIT_CONFIG_GLOBAL": os.devnull,
            "GIT_TERMINAL_PROMPT": "0",
        })
        self.source = Path(__file__).with_name("launch.py").read_text(encoding="utf-8")
        self.git(self.root, "init", "--bare", "--initial-branch=main", str(self.remote))
        self.git(self.root, "init", "--initial-branch=main", str(self.seed))
        (self.seed / "launch.py").write_text(self.source, encoding="utf-8")
        self.write_app(self.seed, "APP_ORIGINAL")
        self.commit(self.seed, "Initial workstation")
        self.git(self.seed, "remote", "add", "origin", str(self.remote))
        self.git(self.seed, "push", "-u", "origin", "main")
        self.git(self.root, "clone", str(self.remote), str(self.checkout))
        self.initial_head = self.git(self.checkout, "rev-parse", "HEAD").stdout.strip()

    def git(self, directory, *args):
        return subprocess.run(
            ["git", "-C", str(directory), *args],
            env=self.env, check=True, text=True, capture_output=True, timeout=15,
        )

    def commit(self, directory, message):
        self.git(directory, "add", ".")
        self.git(directory, "-c", "user.name=Launcher Test", "-c",
                 "user.email=launcher-test@example.invalid", "commit", "-m", message)

    def write_app(self, directory, marker):
        (directory / "app.py").write_text(
            "import json, sys\nprint(" + repr(marker) + ")\n"
            "print('APP_ARGS=' + json.dumps(sys.argv[1:]))\n",
            encoding="utf-8",
        )

    def publish_update(self, update_launcher=False):
        self.write_app(self.seed, "APP_UPDATED")
        if update_launcher:
            (self.seed / "launch.py").write_text(
                self.source.replace("def main():\n", "def main():\n"
                                    "    print('LAUNCHER_UPDATED', flush=True)\n", 1),
                encoding="utf-8",
            )
        self.commit(self.seed, "Publish workstation update")
        self.git(self.seed, "push", "origin", "main")
        return self.git(self.seed, "rev-parse", "HEAD").stdout.strip()

    def launch(self, *args):
        process = subprocess.run(
            [sys.executable, str(self.checkout / "launch.py"), *args],
            cwd=self.checkout, env=self.env, text=True, capture_output=True, timeout=25,
        )
        self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
        return process.stdout

    def head(self):
        return self.git(self.checkout, "rev-parse", "HEAD").stdout.strip()

    def test_clean_checkout_fast_forwards_and_executes_updated_launcher_and_app(self):
        updated_head = self.publish_update(update_launcher=True)
        output = self.launch("--port", "9876")
        self.assertEqual(self.head(), updated_head)
        self.assertIn("LAUNCHER_UPDATED", output)
        self.assertIn("APP_UPDATED", output)
        self.assertIn('APP_ARGS=["--port", "9876"]', output)
        self.assertNotIn("APP_ORIGINAL", output)
        self.assertEqual(output.count("正在检查工作站源码更新"), 1)
        self.assertEqual(self.git(self.checkout, "status", "--porcelain").stdout, "")

    def test_dirty_checkout_preserves_user_files_and_starts_local_app(self):
        self.publish_update()
        self.write_app(self.checkout, "APP_LOCAL_USER_EDIT")
        user_file = self.checkout / "my-sample-groups.csv"
        contents = "SampleID,Group\nA1,CON\n"
        user_file.write_text(contents, encoding="utf-8")
        output = self.launch()
        self.assertEqual(self.head(), self.initial_head)
        self.assertIn("APP_LOCAL_USER_EDIT", output)
        self.assertIn("保留改动并跳过自动更新", output)
        self.assertEqual(user_file.read_text(encoding="utf-8"), contents)

    def test_non_main_branch_and_divergent_main_keep_local_commits(self):
        self.publish_update()
        with self.subTest("non-main branch"):
            self.git(self.checkout, "switch", "-c", "my-analysis")
            self.write_app(self.checkout, "APP_FEATURE_BRANCH")
            self.commit(self.checkout, "Keep custom workflow")
            local_head = self.head()
            output = self.launch()
            self.assertEqual(self.head(), local_head)
            self.assertIn("APP_FEATURE_BRANCH", output)
            self.assertIn("当前不是 main 分支", output)
            self.assertEqual(self.git(self.checkout, "branch", "--show-current").stdout.strip(),
                             "my-analysis")
        with self.subTest("divergent main"):
            self.git(self.checkout, "switch", "main")
            self.write_app(self.checkout, "APP_LOCAL_MAIN_COMMIT")
            self.commit(self.checkout, "Keep local main workflow")
            local_head = self.head()
            output = self.launch()
            self.assertEqual(self.head(), local_head)
            self.assertIn("APP_LOCAL_MAIN_COMMIT", output)
            self.assertIn("自动更新未完成", output)
            self.assertEqual(self.git(self.checkout, "status", "--porcelain").stdout, "")

    def test_unavailable_remote_falls_back_and_no_update_option_skips_check(self):
        self.git(self.checkout, "remote", "set-url", "origin", str(self.root / "missing.git"))
        with self.subTest("unavailable remote"):
            output = self.launch("--port", "9876")
            self.assertEqual(self.head(), self.initial_head)
            self.assertIn("自动更新未完成", output)
            self.assertIn("APP_ORIGINAL", output)
            self.assertIn('APP_ARGS=["--port", "9876"]', output)
        with self.subTest("explicit update opt out"):
            output = self.launch("--no-update", "--port", "9876")
            self.assertEqual(self.head(), self.initial_head)
            self.assertNotIn("正在检查", output)
            self.assertNotIn("自动更新未完成", output)
            self.assertIn('APP_ARGS=["--port", "9876"]', output)


if __name__ == "__main__":
    unittest.main()
