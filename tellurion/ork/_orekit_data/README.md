# Bundled Orekit data

This directory holds a pinned snapshot of
https://gitlab.orekit.org/orekit/orekit-data, populated by
`scripts/update_orekit_data.sh` (the snapshot commit is recorded in
`SNAPSHOT_COMMIT`). If the snapshot is absent, `tellurion.ork.orekit_data`
downloads the data into a user cache at runtime.

With a network connection

```shell
scripts/update_orekit_data.sh
# Alternatively, if the data is in the local cache
# cp -r ~/.cache/tellurion/orekit-data/. tellurion/ork/_orekit_data/
# Check
ls tellurion/ork/_orekit_data | head
git status --short | head
# Update the repository
git add -A tellurion/ork/_orekit_data
git commit -m "Add Orekit data snapshot"
git push

# Rebuild and check
rm -rf dist build *.egg-info
python -m build
tar tzf dist/*.tar.gz | grep tai-utc
unzip -l dist/*.whl | grep tai-utc
```
