test_ID="${1}"
out_fmts=(index number coords features gff bed stats subset table tab subset_table subtab fasta ffn protein faa fa fna contig region)
raw_gff="12383_I10_Bov_IS_ST45.gff"
test_gff="test_${raw_gff}"
if [[ ! -f "${test_gff}" ]] ; then
  cp "${raw_gff}" "${test_gff}"
fi
if [[ -z "${test_ID}" ]] ; then
  test_ID=$(grep -m8 "ID=" "${test_gff}" | tail -1 | cut -d';' -f1 | cut -d'=' -f2)
fi

# Testing basic usage
echo '>>>--' Testing no input '--<<<'
python GF-effortless/proggle.py # should fail with 'no input gff' error
echo -e "\n" '>>>--' Testing no input file '--<<<'
python GF-effortless/proggle.py -f 'ID' -o 'only_exists_if_test_failed.err' # should fail with 'no input gff' error
ls 'only_exists_if_test_failed.err'
echo -e "\n" '>>>--' Testing help menu - see test_helpmenu.txt '--<<<'
python GF-effortless/proggle.py -h > test_helpmenu.txt
echo -e "\n" '>>>--' Testing no parameters '--<<<'
python GF-effortless/proggle.py "${test_gff}" # should reach end of file
echo -e "\n" '>>>--' Testing matchless search '--<<<'
python GF-effortless/proggle.py -s "no matches test." "${test_gff}" # should state no matches found
echo -e "\n" '>>>--' Testing basic search: "${test_ID}" '--<<<'
python GF-effortless/proggle.py -s "${test_ID}" "${test_gff}"

# Testing search functionality
echo -e "\n" '>>>--' Testing search parameters '--<<<'
# 1) basic search
test_gene=$(python GF-effortless/proggle.py -s "${test_ID}" -f "Name" "${test_gff}")
echo "Gene name: - ${test_gene}"
test_family=$(python GF-effortless/proggle.py -s "${test_ID}" -f "Parent" "${test_gff}")
echo "Gene (or Parent) ID: - ${test_family}"
test_family_ID=$(python GF-effortless/proggle.py -s "${test_ID}" -f "locus_tag" "${test_gff}")
echo "locus tag: - ${test_family_ID}"
test_ID=$(python GF-effortless/proggle.py -s "${test_ID}" -f "ID" "${test_gff}")
echo "CDS name: - ${test_ID}"
# 2) search and extract
echo -e "\n" '>>>--' Testing search formatting  - see test_search.txt '--<<<'
for i in ${out_fmts[@]} ; do echo '>>>--' Testing ${i} '--<<<'
  python GF-effortless/proggle.py -s "${test_family_ID}" -f ${i} "${test_gff}"
echo "" ; done > test_search.txt
echo '>>>--' Testing search term "(toxin)" '--<<<' >> test_search.txt
python GF-effortless/proggle.py -s "toxin" -f "Name" "${test_gff}" >> test_search.txt
# 3) search and extract to files
echo -e "\n" '>>>--' Testing search - writing to files '--<<<'
for i in ${out_fmts[@]} ; do outfile=$(echo "${test_gff/.gff}" | sed "s/.gff/_test.${i}/g")
  python GF-effortless/proggle.py -s "${test_family_ID}" -o "${outfile}" -f ${i} "${test_gff}" ; done
echo -e "\n" '>>>--' Testing feature fasta formatting options '--<<<'
echo "Toxins with 60bp upstream:"
python GF-effortless/proggle.py -s "toxin" -f "fasta" -us 60 -o "test_toxins_us.fasta" "${test_gff}"
echo "Toxins with 60AA downstream:"
python GF-effortless/proggle.py -s "toxin" -f "protein" -ds 60 -o "test_toxins_us.aa.fasta" "${test_gff}"
echo "Repeat of above (should raise overwrite error):"
python GF-effortless/proggle.py -s "toxin" -f "protein" -ds 60 -o "test_toxins_us.aa.fasta" "${test_gff}" # should fail without overwrite
echo "Toxins with 60AA downstream, should overwrite successfully with 100AA line splits:"
python GF-effortless/proggle.py -s "toxin" -f "protein" -ds 60 -x -o "test_toxins_us.aa.fasta" -fs 100 "${test_gff}" # should fail without overwrite

# Testing editing functionality
# 1) file-type conversions
echo -e "\n" '>>>--' Testing file type conversions '--<<<'
for i in ${out_fmts[@]} ; do echo '>>>--' Testing conversion to ${i}'--<<<'
python GF-effortless/proggle.py -f ${i} "${test_gff}" ; done
# 2) testing GFF editing
echo -e "\n" '>>>--' Testing file renaming and updating '--<<<'
for i in ${out_fmts[@]} ; do outfile=$(echo "${test_gff/.gff/_r.gff}" | sed "s/.gff/_test.${i}/g")
echo '>>>--' Testing update and conversion to ${i} '--<<<'
python GF-effortless/proggle.py -r -u -o "${outfile}" -f ${i} "${test_gff}" ; done

# confirm no changes to raw GFF_featur
diff -qs "${raw_gff}" "${test_gff}"
