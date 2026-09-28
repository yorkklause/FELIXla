#!/usr/bin/env python3
"""Print the README's upgrade notes for one release tag.

Kept to plain Python so it runs under whatever python3 a runner happens to
have.

The section is found by an explicit `<!-- release-notes: vX.Y.Z -->` anchor
rather than by taking the topmost one, so a release that adds no upgrade notes
falls through to the default rather than silently publishing the previous
release's notes as its own.

Exits 3 when the tag has no section, which is not an error: the caller uses
the default notes. Plain failure is exit 1, which is also what Python uses for
a syntax error or a missing file, so the two must not share a code -- a broken
script that looked like "no section" would publish a plausible release with no
notes at all.
"""

import argparse
import pathlib
import re
import sys

ANCHOR = re.compile(r"^<!--\s*release-notes:\s*(\S+?)\s*-->\s*$")
HEADING = re.compile(r"^#{1,3} ")


def extract(readme, tag):
    lines = readme.split("\n")
    start = None
    for index, line in enumerate(lines):
        match = ANCHOR.match(line)
        if match and match.group(1) == tag:
            start = index + 1
            break
    if start is None:
        return None

    body = []
    for line in lines[start:]:
        # The section runs to the next heading of the same level or higher, or
        # to the next anchor, whichever comes first.
        if ANCHOR.match(line):
            break
        if HEADING.match(line) and body:
            break
        body.append(line)

    # Drop the section's own heading; the release page supplies the title.
    while body and (not body[0].strip() or HEADING.match(body[0])):
        body.pop(0)
    while body and not body[-1].strip():
        body.pop()
    return "\n".join(body) if body else None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("tag")
    parser.add_argument("readme", type=pathlib.Path)
    args = parser.parse_args()

    notes = extract(args.readme.read_text(), args.tag)
    if notes is None:
        print("no release notes section for %s" % args.tag, file=sys.stderr)
        return 3
    print(notes)
    return 0


if __name__ == "__main__":
    sys.exit(main())
