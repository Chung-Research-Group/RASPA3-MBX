"""Exercise the release workflow against local repositories without network access."""

import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SYNC = Path(__file__).resolve().parents[1] / ".github/scripts/sync_release.py"


class ReleaseSyncTest(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.source = self.root / "upstream"
        self.worker = self.root / "worker"
        self.downstream = self.root / "downstream.git"
        self.releases = self.root / "stable"
        self.output = self.root / "outputs"
        self.env = {
            **os.environ,
            "GIT_AUTHOR_NAME": "RASPA test",
            "GIT_COMMITTER_NAME": "RASPA test",
            "GIT_AUTHOR_EMAIL": "test@example.invalid",
            "GIT_COMMITTER_EMAIL": "test@example.invalid",
            "GITHUB_OUTPUT": str(self.output),
        }
        self.git(self.root, "init", "-b", "main", str(self.source))
        self.commit(self.source, "shared", "original")
        self.git(self.source, "tag", "v3.0.30")
        self.git(self.root, "clone", "--bare", str(self.source), str(self.downstream))
        self.git(self.root, "clone", str(self.downstream), str(self.worker))
        self.git(self.worker, "push", "origin", "HEAD:refs/heads/mbx")
        self.releases.write_text("2026-05-31T07:44:20Z\tv3.0.30\n")

    def git(self, cwd, *args):
        return subprocess.run(
            ["git", *args], cwd=cwd, env=self.env, check=True, capture_output=True, text=True
        ).stdout.strip()

    def commit(self, cwd, file, content):
        (cwd / file).write_text(content)
        self.git(cwd, "add", file)
        self.git(cwd, "commit", "-m", content)
        return self.git(cwd, "rev-parse", "HEAD")

    def run_sync(self, *args, success=True):
        result = subprocess.run(
            [sys.executable, str(SYNC), *args], cwd=self.worker,
            env=self.env, capture_output=True, text=True
        )
        self.assertEqual(result.returncode == 0, success, result.stderr)

    def prepare(self, success=True):
        self.output.write_text("")
        self.run_sync("prepare", str(self.releases), "--upstream", str(self.source), success=success)
        return dict(line.split("=", 1) for line in self.output.read_text().splitlines())

    def promote(self, candidate, success=True):
        self.run_sync("promote", *(candidate[k] for k in ("base", "branch", "sha", "release")),
                      success=success)

    def release(self, tag, path, content):
        sha = self.commit(self.source, path, content)
        self.git(self.source, "tag", "-a", tag, "-m", tag)
        with self.releases.open("a") as stream:
            stream.write(f"2026-11-19T01:00:00Z\t{tag}\n")
        return sha

    def tip(self, branch):
        return self.git(self.downstream, "rev-parse", "refs/heads/" + branch)

    def test_only_published_tags_are_merged_in_order(self):
        self.assertEqual(self.prepare(), {"sha": ""})
        first = self.release("v3.0.31", "first", "release one")
        second = self.release("v3.0.32", "second", "release two")
        self.commit(self.source, "unreleased", "later upstream work")
        candidate = self.prepare()
        self.assertEqual(candidate["release"], first)
        self.assertEqual(self.tip("main"), candidate["base"])
        self.assertEqual(self.prepare(), candidate)
        self.promote(candidate)
        self.assertEqual(self.tip("main"), self.tip("mbx"))
        self.assertEqual(self.tip("main"), self.tip(candidate["branch"]))
        candidate = self.prepare()
        self.assertEqual(candidate["release"], second)
        self.promote(candidate)
        self.assertEqual(self.prepare(), {"sha": ""})
        tree = self.git(self.downstream, "ls-tree", "--name-only", "main").splitlines()
        self.assertNotIn("unreleased", tree)

    def test_dry_run_leaves_all_remote_refs_untouched(self):
        self.release("v3.0.31", "new", "new release")
        refs = self.git(self.downstream, "show-ref")
        self.run_sync("prepare", str(self.releases), "--dry-run", "--upstream", str(self.source))
        self.assertEqual(self.git(self.downstream, "show-ref"), refs)
        self.assertIn("tag=v3.0.31", self.output.read_text())

    def test_conflict_and_manual_resolution(self):
        self.commit(self.worker, "shared", "fork revision")
        self.git(self.worker, "push", "origin", "HEAD:main")
        upstream = self.release("v3.0.31", "shared", "upstream revision")
        self.prepare(success=False)
        self.assertEqual(self.tip("sync-to-upstream/v3.0.31"), self.tip("main"))
        self.git(self.worker, "fetch", str(self.source), "refs/tags/v3.0.31")
        self.git(self.worker, "merge", "--no-edit", "-X", "ours", upstream)
        self.git(self.worker, "push", "origin", "HEAD:refs/heads/sync-to-upstream/v3.0.31")
        self.promote(self.prepare())
        self.assertEqual(self.prepare(), {"sha": ""})

    def test_head_change_requires_fresh_validation(self):
        self.release("v3.0.31", "new", "release")
        original = self.prepare()
        advanced = self.commit(self.worker, "fork", "downstream work")
        self.git(self.worker, "push", "origin", "HEAD:main")
        self.promote(original, success=False)
        self.assertEqual(self.tip("main"), advanced)
        updated = self.prepare()
        self.assertNotEqual(original["sha"], updated["sha"])
        self.promote(updated)
        self.assertEqual(self.tip("main"), self.tip("mbx"))

    def test_diverged_mbx_is_preserved(self):
        self.release("v3.0.31", "new", "release")
        original = self.prepare()
        other = self.commit(self.worker, "fork", "independent mbx work")
        self.git(self.worker, "push", "origin", "HEAD:mbx")
        self.promote(original, success=False)
        self.assertEqual(self.tip("mbx"), other)
        self.assertEqual(self.tip("main"), original["base"])


if __name__ == "__main__":
    unittest.main()
