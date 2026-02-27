# Frequently Asked Questions (FAQ)

## My search is taking a long time. What can I do?

Most likely your search is returning a lot of results. The search commands have
several options to reduce the number of results returned, such as `--limit`.

## How can I run workflows with many structure files on an HPC filesystem?

If your cluster home/project filesystem is slow with many small files, run the
high-file-count steps in local scratch and archive results to persistent
storage. A common scratch policy is automatic cleanup of files older than 14
days, but this policy is cluster-specific.

Use the same commands from the README, but point temporary directories and cache
to scratch:

```shell
SCRATCH_RUN="${SCRATCH:-/tmp}/${USER}/protein-quest-run"
PERSIST_DIR="$HOME/protein-quest-results"
mkdir -p "$SCRATCH_RUN/cache" "$PERSIST_DIR"

# Example: retrieve and filter on scratch (many files)
protein-quest retrieve alphafold \
	--cache-dir "$SCRATCH_RUN/cache" \
	alphafold.csv "$SCRATCH_RUN/downloads-af"

protein-quest filter confidence \
	--confidence-threshold 50 \
	--min-residues 100 \
	--max-residues 1000 \
	--cache-dir "$SCRATCH_RUN/cache" \
	"$SCRATCH_RUN/downloads-af" "$SCRATCH_RUN/filtered"

# Bundle outputs for persistent storage (home/project)
tar -C "$SCRATCH_RUN" -cf "$PERSIST_DIR/filtered-$(date +%F).tar" filtered
```

To resume later, unpack the tarball again to scratch:

```shell
mkdir -p "$SCRATCH_RUN"
tar -C "$SCRATCH_RUN" -xf "$PERSIST_DIR/filtered-2026-02-27.tar"
```

This pattern keeps active runs fast on scratch while storing long-term results
as a small number of tarballs in persistent storage.

If you have downstream tools that can not read from tarballs, you can mount the
tarball as a filesystem with [ratarmount](https://github.com/mxmlnkn/ratarmount)
and unmount with `fusermount -u <mount_point>`.

## My log is polluted with progress bar lines. How can I fix this?

To reduce the number of lines printed by the progress bar, you can increase the
minimum interval between updates with the `TQDM_MININTERVAL` environment
variable. For example, setting it to `9` will update the progress bar every 9
seconds instead of every 0.1 seconds.

To not have any progress bars at all, you can set `TQDM_DISABLE` environment
variable to any value.

## My protein-quest question is not answered here. Where can I get help?

Please see the [Contributing](CONTRIBUTING.md#you-have-a-question) document for
instructions on how to ask questions and report issues.
