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

If you have downstream tools that cannot read from tarballs, you can mount the
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

## I do not like the default values for some parameters. How can I change them?

The simple approach is to specify the parameter values on the command line. If
you do not want to specify them every time, you can use a configuration file to
override the default values.

The user configuration file is located at `~/.config/protein-quest/config.toml`.
See
[determine_config_file_location](https://www.bonvinlab.org/protein-quest/autoapi/protein_quest/cli/common.md#protein_quest.cli.common.determine_config_file_location)
for more places where the configuration file can be located.

For example, to set the default value for the `--confidence-threshold` parameter
of the `filter confidence` command to 50, you can add the following to your
`config.toml` file:

```toml
[filter.confidence]
confidence_threshold = 50
```

The default values displayed in `--help` will remain unchanged, but the values
in the configuration file will be used when you run commands and do not specify
those parameters on the command line.

> Configuration file does not play well with RO Crate provenance recording
>
> A generated RO Crate (using `--prov` parameter) only stores the parameters
> used on the command line. It does not keep track of the config file or its
> contents. Therefore, using a config file for parameters that affect the
> reproducibility of your workflow is not recommended.

## On the same machine multiple users want to share the same cache. How can we do that?

You can set a shared cache directory that all users have read and write access
to.

1. First make sure all users are member of the same group, for example `pq`.
2. Create a shared cache directory and set the group ownership and permissions
   appropriately. For example with
   `mkdir -p /shared/protein-quest-cache;chgrp -R pq /shared/protein-quest-cache;chmod -R g+rwxs,o-rwx /shared/protein-quest-cache`.
3. Create a central configuration file located at
   `/shared/protein-quest-cache/config.toml` with content like
   [example-shared-config.toml](example-shared-config.toml) file.
4. Ask users to put in their shell startup file the following line:
   `PROTEIN_QUEST_CONFIG=/shared/protein-quest-cache/config.toml` to ensure that
   shared cache is used. Alternatively, users can symlink the shared config file
   to their home directory as `~/.config/protein-quest/config.toml`.
5. Also ask users to put in their shell startup file the following line:
   `umask u=,g=,o=rwx` to ensure that new files can be read/written by yourself
   and group members (and no access to others).
