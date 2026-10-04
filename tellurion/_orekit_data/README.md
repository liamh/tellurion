# Bundled Orekit data

This directory holds a pinned snapshot of
https://gitlab.orekit.org/orekit/orekit-data, populated by
`scripts/update_orekit_data.sh` (the snapshot commit is recorded in
`SNAPSHOT_COMMIT`). If the snapshot is absent, `tellurion.orekit_data`
downloads the data into a user cache at runtime.
