"""Prepare and promote tested downstream branches for an official RASPA3 release."""

import argparse
import os
from pathlib import Path
import subprocess
import sys


def git(*args):
    result = subprocess.run(["git", *args], text=True, stdout=subprocess.PIPE, check=True)
    return result.stdout.strip()


def ancestor(older, newer):
    status = subprocess.run(["git", "merge-base", "--is-ancestor", older, newer]).returncode
    if status not in (0, 1):
        raise RuntimeError(f"Cannot compare {older} and {newer}")
    return status == 0


def remote_tip(branch):
    ref = "refs/heads/" + branch
    line = git("ls-remote", "--refs", "origin", ref)
    if not line:
        return ""
    git("fetch", "--no-tags", "origin", ref)
    return git("rev-parse", "FETCH_HEAD^{commit}")


def merge(base, other, message):
    if ancestor(other, base):
        return base
    tree = git("merge-tree", "--write-tree", base, other)
    return git("commit-tree", tree, "-p", base, "-p", other, "-m", message)


def outputs(**values):
    with Path(os.environ["GITHUB_OUTPUT"]).open("a") as stream:
        for key, value in values.items():
            print(f"{key}={value}", file=stream)
            print(f"{key}={value}")


def prepare(args):
    base = remote_tip("main")
    if not base:
        raise RuntimeError("Downstream main does not exist")
    # The backup includes v3.0.30. Older publications predate this workflow.
    releases = sorted(line.split("\t") for line in args.releases.read_text().splitlines())
    for date, tag in releases:
        if date < "2026-05-31T07:44:20Z":
            continue
        ref = "refs/tags/" + tag
        git("check-ref-format", ref)
        git("fetch", "--no-tags", args.upstream, ref)
        release = git("rev-parse", "FETCH_HEAD^{commit}")
        if ancestor(release, base):
            continue

        branch = "sync-to-upstream/" + tag
        candidate = remote_tip(branch)
        if not candidate:
            # Publish the starting point first so conflicts have a branch for repair.
            if not args.dry_run:
                git("push", "--force-with-lease=refs/heads/" + branch + ":",
                    "origin", base + ":refs/heads/" + branch)
            candidate = base
        try:
            tested = merge(candidate, base, f"Merge current downstream main into {branch}")
            tested = merge(tested, release, f"Merge official RASPA3 {tag} into downstream")
        except subprocess.CalledProcessError as error:
            conflicts = (error.stdout or "").strip()
            raise RuntimeError(f"Resolve the merge on {branch} and rerun.\n{conflicts}") from error

        if tested != candidate and not args.dry_run:
            git("push", "--force-with-lease=refs/heads/" + branch + ":" + candidate,
                "origin", tested + ":refs/heads/" + branch)
        outputs(sha=tested, base=base, branch=branch, release=release, tag=tag)
        return
    print("All published upstream releases are present in downstream main.")
    outputs(sha="")


def promote(args):
    if remote_tip(args.branch) != args.sha or remote_tip("main") != args.base:
        raise RuntimeError("Downstream main or the sync branch changed. Rerun for a fresh validation.")
    if not ancestor(args.base, args.sha) or not ancestor(args.release, args.sha):
        raise RuntimeError("Validated commit must contain the downstream base and the upstream release")
    # If identified as a release alias, mbx must be advanced with main, never overwritten.
    alias = remote_tip("mbx")
    leases = ["--force-with-lease=refs/heads/main:" + args.base]
    refs = [args.sha + ":refs/heads/main"]
    if alias:
        if not ancestor(alias, args.sha):
            raise RuntimeError("Downstream mbx diverged. Reconcile the two combined branches first.")
        leases.append("--force-with-lease=refs/heads/mbx:" + alias)
        refs.append(args.sha + ":refs/heads/mbx")
    git("push", "--atomic", *leases, "origin", *refs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    stage = sub.add_parser("prepare")
    stage.add_argument("releases", type=Path)
    stage.add_argument("--dry-run", action="store_true")
    stage.add_argument("--upstream", default="https://github.com/iRASPA/RASPA3.git")
    finish = sub.add_parser("promote")
    for name in ("base", "branch", "sha", "release"):
        finish.add_argument(name)
    args = parser.parse_args()
    if args.action == "prepare":
        prepare(args)
    else:
        promote(args)


if __name__ == "__main__":
    try:
        main()
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        sys.exit(str(error))
