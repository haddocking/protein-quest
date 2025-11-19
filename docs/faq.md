# Frequently Asked Questions (FAQ)

## My search is taking a long time. What can I do?

Most likely your search is returning a lot of results.
The search commands have several options to reduce the number of results returned, such as `--limit`.

## My log is polluted with progress bar lines. How can I fix this?

To reduce the number of lines printed by the progress bar, you can increase the minimum interval between updates with the `TQDM_MININTERVAL` environment variable.
For example, setting it to `9` will update the progress bar every 9 seconds instead of every 0.1 seconds.

To not have any progress bars at all, you can set `TQDM_DISABLE` environment variable to any value.

## My protein-quest question is not answered here. Where can I get help?

Please see the [Contributing](CONTRIBUTING.md#you-have-a-question) document for instructions on how to ask questions and report issues.