// Version string for FELIXla, injected by the Makefile from `git describe` so
// a prefix records the exact tree it came from, including `-dirty` for a build
// made from uncommitted work. Source tarballs with no git metadata fall back to
// the literal below, which must be bumped alongside a release tag.
#ifndef FELIXLA_VERSION_H
#define FELIXLA_VERSION_H

#ifndef FELIXLA_VERSION
#define FELIXLA_VERSION "0.5.0-unknown"
#endif

#endif
