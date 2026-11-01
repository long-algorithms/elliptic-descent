# Release publication and draft recovery

## Draft-release recovery

The publisher lists authenticated releases and retains the numeric release ID.
Draft creation, metadata updates and asset uploads use that ID rather than a
release-by-tag lookup. Every metadata write explicitly supplies both the intended
`tag_name` and `target_commitish`, as well as publication state. An omitted
`tag_name` during the first recovery changed the draft to a GitHub-generated
`untagged-...` placeholder; a stateful regression reproduces that failure.

Before publication, a fully paginated asset listing is checked against every
freshly audited file's size and SHA256. Identical draft assets may be reused;
changed expected files are replaced by numeric asset ID. Unknown or duplicate
assets fail instead of being silently deleted. Missing, changed or incomplete
uploads prevent publication.

An abandoned draft targeting an older revision is not replaced implicitly.
The narrowly scoped `--replace-draft ID@OLD_SHA` option requires the exact
draft ID, an unpublished/non-prerelease state, no existing tag, and proof that
OLD_SHA is an ancestor of the newly audited revision. The full exact-main
workflow and artifact audits must pass before any draft mutation. Published
releases and tags are never retargeted. The release workflow records the one
known draft left by the initial publication failure; the allowance has no
effect once the same-revision release has been published.

A draft whose intended tag was replaced by a placeholder additionally requires
`--recover-draft-tag EXACT_OLD_TAG`. Recovery identifies it by the supplied
numeric ID and checks the exact old revision, exact old tag, unpublished state,
absence of both Git tags, and ancestry. No wildcard `untagged-*` match is used.
The workflow records the observed triple for draft 430578096. All exact-main
mathematical, package and evidence gates still run before any mutation. A retry
after successful retargeting uses the ordinary intended-tag draft path.

## Every exact-main gate must run

All nine required audit workflows run on every push to `main`,
including publication-tooling or documentation-only changes, and on every pull
request. They no longer rely on source-path changes to satisfy the exact-main
publication inventory. Duplicate feature-branch push runs are avoided; every
pull-request revision still receives the full required matrix.

Group construction and portable replay remain separate mandatory computations,
now in dependent jobs with a 20-minute budget each. A runner-dependent combined
timeout cannot be treated as a verified replay. The construction job uploads
`.record` transfer files; only the replay job produces the final fixture and
replay filenames consumed by the unchanged release gate. Every one of the eight
fixtures must pass both stages. This is job isolation, not a reduced verifier,
a cached mathematical verdict, or an increased per-job timeout.
