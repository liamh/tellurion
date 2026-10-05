# Release Checklist for Tellurion

## Pre-Release (on `develop`)

- [ ] All feature branches merged and tests passing
- [ ] Update `version` in `pyproject.toml` (e.g., `0.1.0` → `0.2.0`)
- [ ] Update `CHANGELOG.md` or release notes with new features, fixes, and breaking changes
- [ ] Commit these changes to `develop` with message: `"Bump version to X.Y.Z"`
- [ ] Create a pull request from `develop` to `production` for review

## Release (on `production`)

- [ ] Merge the PR to `production`
- [ ] Verify the merge commit is on `production`: `git log origin/production --oneline -1`
- [ ] Create an annotated tag: `git tag -a vX.Y.Z -m "Release X.Y.Z"` (e.g., `v0.2.0`)
- [ ] Push the tag: `git push origin vX.Y.Z`
- [ ] Monitor the GitHub Actions tab: the `publish.yaml` workflow should trigger
- [ ] Verify the workflow completes successfully (check logs if it fails)
- [ ] Confirm the new version appears on [PyPI](https://pypi.org/project/tellurion/)

## Post-Release (on `develop`)

- [ ] Verify PyPI installation works: `pip install --upgrade tellurion`
- [ ] Create a GitHub Release from the tag on GitHub (optional but recommended for visibility)
- [ ] Announce the release if applicable

## Docs-Only Update (no release)

Documentation or metadata changes that don't affect code can be applied to `production` without tagging a release.

```bash
git checkout production
git merge develop        # or cherry-pick specific commits
git push origin production
```

**Important notes:**
- The workflow only runs on `v*` tags, so no PyPI publish occurs.
- If the change is in `README.md` (your PyPI package description), PyPI will not update until the next release.
- If you need PyPI to reflect a docs change, publish a post-release: bump `version` to `X.Y.Z.post1` in `pyproject.toml`, then tag and push as normal.
- If using cherry-pick to bring only specific commits, be careful not to accidentally mix in unreleased code changes.
- Check your Read the Docs configuration to see which branch it builds from (may be `production`, `develop`, or tags only).

## Troubleshooting

**Workflow did not trigger**
- Verify the tag was pushed: `git push origin vX.Y.Z`
- Check that the tag points to a commit on `production`: `git log vX.Y.Z --oneline -1` and verify ancestry with `git merge-base --is-ancestor vX.Y.Z origin/production`

**Workflow failed at "Require tagged commit to be on production"**
- The tag is not reachable from `production`. Ensure you tagged a commit that was merged to `production`.

**PyPI publish step failed**
- Verify the PyPI trusted publisher is configured with:
  - Owner: `liamh`
  - Repository: `tellurion`
  - Workflow: `publish.yaml`
  - Environment: `pypi`
- Check the workflow logs for the full error message.

**Need to re-release (e.g., tag was pushed to wrong commit)**
- Delete the tag locally and remotely: `git tag -d vX.Y.Z` and `git push origin --delete vX.Y.Z`
- Fix the issue (re-merge, update version, etc.) and try again
