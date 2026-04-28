#!/bin/bash
# File: audit_pdb2cif_fragility.sh · Git: Michael G. Prisant (uncommitted) · Txd: 04/28/26 04:13 with local modifications!
#
# Read-only audit: does `pdb2cif | mkdssp --output-format=dssp` handle every
# fragility category that pydangle's _clean_pdb_for_dssp currently patches?
#
# Tests one or two representative PDBs from each fragility class:
#   BLANK-SEQRES-COL12       — blank chain ID in SEQRES col 12 (top100)
#   MOD-AA-HETATM            — modified amino acids deposited as HETATM (top100)
#   HETATM-SEQRES-COLLISION  — HOH/ligand HETATM colliding with SEQRES (top500)
#   HEADER-COMPND-OVERFLOW   — HEADER/COMPND fields exceed legacy width (top8000)
#   CHAIN-OVERFLOW-GT26      — input has >26 case-sensitive chain IDs (top2018)
#
# Per case: input PDB chain set vs output classic-DSSP chain set. PASS = no chain
# lost between input and output. Reported alongside the fragility markers found
# in the input (ANISOU/SIGATM/SIGUIJ, blank chain col 12/22, HETATM, NCHAINS).
#
# Result of 2026-04-28 run: 4 PASS / 12 FAIL of 16 cases. Reported in
# session_summary_2026-04-28.md. Confirms pdb2cif is NOT a panacea —
# (a) cleanly handles the modified-AA HETATM class but (b) introduces 5 new
# strict-validation failures (blank-chain SEQRES, pre-Y2K dates, duplicate
# REVDAT keys) and (c) does not solve the chain-overflow class because the
# 26-chain cap lives in mkdssp's classic-DSSP output writer, not its PDB-input
# parser.
#
# Usage: ./audit_pdb2cif_fragility.sh
# Inputs: /home/prisant/Desktop/Data/{top100_pdbs,top500_pdbs,top8000_pdbs_hom70,top2018_pdbs_*_filtered_hom70}
# Workdir: /tmp/pdb2cif_audit (cif/dssp/logs preserved for inspection)

WORK=/tmp/pdb2cif_audit
rm -rf "$WORK" && mkdir -p "$WORK"

# (label, dataset, pdbid, fragility-category)
CASES=(
  "BLANK-SEQRES-COL12        | top100_pdbs                       | 3b5c | blank chain in SEQRES col 12"
  "BLANK-SEQRES-COL12        | top100_pdbs                       | 4ptp | blank chain in SEQRES col 12"
  "MOD-AA-HETATM             | top100_pdbs                       | 1cpc | PCA/MEN as HETATM (must keep)"
  "MOD-AA-HETATM             | top100_pdbs                       | 1lkk | PTR as HETATM (must keep)"
  "MOD-AA-HETATM             | top100_pdbs                       | 1smd | PCA as HETATM (must keep)"
  "MOD-AA-HETATM             | top100_pdbs                       | 2er7 | modified-AA HETATM (must keep)"
  "HETATM-SEQRES-COLLISION   | top500_pdbs                       | 1a1y | HOH/ligand collision (silent null)"
  "HETATM-SEQRES-COLLISION   | top500_pdbs                       | 1gdo | HOH/ligand collision (silent null)"
  "HETATM-SEQRES-COLLISION   | top500_pdbs                       | 1tax | HOH/ligand collision (silent null)"
  "HEADER-COMPND-OVERFLOW    | top8000_pdbs_hom70                | 1h64 | header field overflow"
  "HEADER-COMPND-OVERFLOW    | top8000_pdbs_hom70                | 1ryp | header field overflow"
  "HEADER-COMPND-OVERFLOW    | top8000_pdbs_hom70                | 2zzs | header field overflow"
  "HEADER-COMPND-OVERFLOW    | top8000_pdbs_hom70                | 3mt6 | header field overflow"
  "CHAIN-OVERFLOW-GT26       | top2018_pdbs_full_filtered_hom70  | 3uoi | 48 chains"
  "CHAIN-OVERFLOW-GT26       | top2018_pdbs_mc_filtered_hom70    | 4ub8 | 39 chains"
  "CHAIN-OVERFLOW-GT26       | top2018_pdbs_mc_filtered_hom70    | 5b66 | 38 chains"
)

# Returns the set of input ATOM chain IDs (case-sensitive, single char), one per line.
input_chains() {
  awk 'substr($0,1,6)=="ATOM  " && length($0)>=22 {ch=substr($0,22,1); if (ch!=" ") print ch}' "$1" | sort -u
}

# Returns the set of output DSSP chain IDs (column 12, case-sensitive).
output_chains() {
  awk '
    /^  #  RESIDUE/ {f=1; next}
    f && length($0)>=14 {
      aa=substr($0,14,1); if (aa=="!") next
      ch=substr($0,12,1); if (ch!=" ") print ch
    }
  ' "$1" | sort -u
}

# Detect fragility markers actually present in the input.
fragility_markers() {
  awk '
    BEGIN { anisou=0; sigatm=0; siguij=0; blank22=0; blank12=0; hetatm=0; chains[""]=0 }
    /^ANISOU/  { anisou=1 }
    /^SIGATM/  { sigatm=1 }
    /^SIGUIJ/  { siguij=1 }
    /^HETATM/  { hetatm=1 }
    /^ATOM  / && length($0)>=22 {
      ch=substr($0,22,1)
      if (ch==" ") blank22=1
      chains[ch]=1
    }
    /^SEQRES/ && length($0)>=12 {
      ch=substr($0,12,1)
      if (ch==" ") blank12=1
    }
    END {
      n=0; for (c in chains) if (c!="") n++
      tags=""
      if (anisou)  tags=tags"ANISOU "
      if (sigatm)  tags=tags"SIGATM "
      if (siguij)  tags=tags"SIGUIJ "
      if (hetatm)  tags=tags"HETATM "
      if (blank22) tags=tags"BLANK-COL22 "
      if (blank12) tags=tags"BLANK-COL12 "
      if (n>26)    tags=tags"NCHAINS="n" "
      else         tags=tags"nchains="n" "
      print tags
    }
  ' "$1"
}

printf "%-26s | %-6s | %-6s | %-9s | %-9s | %-7s | %s\n" \
       "Category" "PDBID" "Result" "InChains" "OutChains" "Lost" "Markers / Notes"
printf "%s\n" "------------------------------------------------------------------------------------------------------------------------------"

PASS=0; FAIL=0; TOTAL=0
for entry in "${CASES[@]}"; do
  TOTAL=$((TOTAL+1))
  IFS='|' read -r label ds pid note <<<"$entry"
  label=$(echo "$label" | xargs); ds=$(echo "$ds" | xargs); pid=$(echo "$pid" | xargs); note=$(echo "$note" | xargs)
  pre="${pid:0:2}"
  raw="/home/prisant/Desktop/Data/$ds/$pre/$pid/$pid.pdb"
  cif="$WORK/$pid.cif"; dssp="$WORK/$pid.dssp"; logf="$WORK/$pid.log"

  if [ ! -f "$raw" ]; then
    printf "%-26s | %-6s | MISSING-INPUT\n" "$label" "$pid"; FAIL=$((FAIL+1)); continue
  fi

  in_set=$(input_chains "$raw"); n_in=$(echo "$in_set" | grep -c .)
  markers=$(fragility_markers "$raw")

  pdb2cif "$raw" "$cif" >"$logf" 2>&1
  rc_p2c=$?
  if [ $rc_p2c -ne 0 ] || [ ! -s "$cif" ]; then
    printf "%-26s | %-6s | %-6s | %-9s | %-9s | %-7s | pdb2cif rc=%d\n" \
           "$label" "$pid" "FAIL" "$n_in" "-" "-" "$rc_p2c"
    FAIL=$((FAIL+1)); continue
  fi

  mkdssp --output-format=dssp "$cif" "$dssp" >>"$logf" 2>&1
  rc_dssp=$?
  if [ ! -s "$dssp" ]; then
    printf "%-26s | %-6s | %-6s | %-9s | %-9s | %-7s | mkdssp empty rc=%d\n" \
           "$label" "$pid" "FAIL" "$n_in" "-" "-" "$rc_dssp"
    FAIL=$((FAIL+1)); continue
  fi

  out_set=$(output_chains "$dssp"); n_out=$(echo "$out_set" | grep -c .)
  lost=$(comm -23 <(echo "$in_set") <(echo "$out_set") | tr '\n' ',' | sed 's/,$//')
  n_lost=$(echo "$lost" | awk -F, 'NF{print NF; exit}')
  [ -z "$n_lost" ] && n_lost=0

  if [ "$n_lost" = "0" ]; then
    res="PASS"; PASS=$((PASS+1))
  else
    res="FAIL"; FAIL=$((FAIL+1))
  fi

  printf "%-26s | %-6s | %-6s | %-9s | %-9s | %-7s | %s\n" \
         "$label" "$pid" "$res" "$n_in" "$n_out" "${lost:-none}" "$markers"
done

echo
echo "Summary: $PASS PASS / $FAIL FAIL  (of $TOTAL)"
echo "Workdir: $WORK   (cif/dssp/logs preserved for any failure inspection)"
