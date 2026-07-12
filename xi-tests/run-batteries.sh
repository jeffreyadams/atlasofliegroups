#!/bin/bash
# Regenerate the xi-tests stored outputs. Run from the atlasCovering root.
#
# HARNESS NOTE: the cover batteries use the LIGHT preload below (basic +
# hermitian + induction). Do NOT preload all.at for these: since all.at
# gained the FPP suite (to_ht overrides of full_deform etc.), the FPP layer
# breaks deformation/hermitian calls on COVER parameters (its hash/time
# machinery and delta-fixed assumptions do not know about covers). That
# interaction is documented in coverNotes/e2-covered-endoscopy.tex and is a
# known open item; plain-group (xi=0) behavior under all.at is unaffected.
PRELOAD="atlas-scripts/basic.at atlas-scripts/hermitian.at atlas-scripts/induction.at"
for t in cover-smoke gate2 hermitian-cover tensor-bijection induction-cover; do
  [ -f xi-tests/$t.in ] || continue
  ./atlas --path=atlas-scripts $PRELOAD < xi-tests/$t.in 2>&1 | sed -n '/^###/,$p' > xi-tests/$t-out.txt
  echo "regenerated $t-out.txt ($(wc -l < xi-tests/$t-out.txt) lines)"
done
# endoscopy batteries have their own preload
./atlas --path=atlas-scripts atlas-scripts/endoscopy_xi.at atlas-scripts/groups.at < xi-tests/endoscopy-xi.in 2>&1 | sed -n '/^###/,$p' > xi-tests/endoscopy-xi-out.txt
echo "regenerated endoscopy-xi-out.txt"
./atlas --path=atlas-scripts atlas-scripts/endoscopy_xi.at atlas-scripts/groups.at < xi-tests/endoscopy-lift-xi.in 2>&1 | sed -n '/^###/,$p' > xi-tests/endoscopy-lift-xi-out.txt
echo "regenerated endoscopy-lift-xi-out.txt"
