# Release publication and draft recovery

## Draft-release recovery

The publisher lists authenticated releases and retains the numeric release ID;
it does not assume an unpublished release is available from the tag endpoint.
Before publication it compares every uploaded asset's size and SHA256 digest
with the freshly audited local artifact. Duplicate, missing or changed assets
prevent publication.

An abandoned draft targeting an older revision is not replaced implicitly.
The narrowly scoped `--replace-draft ID@OLD_SHA` option requires the exact
draft ID, an unpublished/non-prerelease state, no existing tag, and proof that
OLD_SHA is an ancestor of the newly audited revision. The full exact-main
workflow and artifact audits must pass before any draft mutation. Published
releases and tags are never retargeted. The release workflow records the one
known draft left by the initial publication failure; the allowance has no
effect once the same-revision release has been published.
