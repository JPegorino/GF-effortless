#!/bin/bash
# This script was created for author testing and may not work on your system.
# A single GFF file, "UnitTestGenome.gff" is provided for basic testing purposes.

# Read User Inputs
while [ "$1" != "" ]; do
  case $1 in
    -g | --gff )	shift
			raw_gff=$1
      ;;
	  -i | --locus_ID )		shift
			test_ID=$1
			;;
		-h | --help )
          printf "\nUsage: Use unit_test.sh [-g GFF_FILEPATH] [-i LOCUS_ID]\n"
    			exit 0
      ;;
    * )                     printf "\nUnrecognised option:\tUse unit_test.sh -h | unit_test.sh --help\n"
    exit 1
  esac
  shift
done

# Parse User Inputs
if [[ ! -f "${raw_gff}" ]] ; then
  raw_gff="UnitTestGenome.gff"
fi
test_gff="test_$(basename ${raw_gff})"
if [[ ! -f "${test_gff}" ]] ; then
  cp "${raw_gff}" "${test_gff}" || exit 1
fi
if [[ -z "${test_ID}" ]] ; then
  test_ID=$(grep -m8 "ID=" "${test_gff}" | tail -1 | cut -d';' -f1 | cut -d'=' -f2)
fi

# Prepare Software & Diagnostics Logs
out_fmts=(index number coords features gff bed stats subset table tab subset_table subtab fasta ffn protein faa fa fna)
diagnostics_log="test_run_unittest.log"
if [[ -f "${diagnostics_log}" ]] ; then
  echo 'Diagnostics log already exists... overwrite?'
  rm -i "${diagnostics_log}"
fi
if [[ -f "${diagnostics_log}" ]] ; then
  echo "Okay, will not overwrite. Exiting..."
  exit 1
fi
echo "Running Proggle diagnostics..."
exec 1> "$diagnostics_log" 2>&1

# Testing basic usage - exit if any of these do not succeed but should
echo '>>>--' Testing help menu - see test_helpmenu.txt '--<<<'
python ../proggle.py -h > test_helpmenu.txt || exit 1 
echo -e "\n" '>>>--' Testing no input '--<<<'
python ../proggle.py # should fail with 'no input gff' error
echo -e "\n" '>>>--' Testing no input file '--<<<'
python ../proggle.py -f 'ID' -o 'only_exists_if_test_failed.err' # should fail with 'no input gff' error
ls 'only_exists_if_test_failed.err'
echo -e "\n" '>>>--' Testing no parameters '--<<<'
python ../proggle.py "${test_gff}" || exit 1  # should reach end of file
echo -e "\n" '>>>--' Testing basic search: "${test_ID}" '--<<<'
python ../proggle.py -s "${test_ID}" "${test_gff}"
echo -e "\n" '>>>--' Testing matchless search '--<<<'
python ../proggle.py -s "no matches test." "${test_gff}" # should state no matches found

# Testing search functionality
echo -e "\n\n" '>>>--' Testing search parameters '--<<<'
###########
# 1) basic search
test_gene=$(python ../proggle.py -s "${test_ID}" -f "Name" "${test_gff}")
echo "Gene name: - ${test_gene}"
test_family=$(python ../proggle.py -s "${test_ID}" -f "Parent" "${test_gff}")
echo "Gene (or Parent) ID: - ${test_family}"
test_family_ID=$(python ../proggle.py -s "${test_ID}" -f "locus_tag" "${test_gff}")
echo "locus tag: - ${test_family_ID}"
test_ID=$(python ../proggle.py -s "${test_ID}" -f "ID" "${test_gff}")
echo "CDS name: - ${test_ID}"
# 2) search and extract
out_fmts=(index number coords bed stats subset table tab subset_table subtab fasta ffn protein faa fa fna)
echo -e "\n" '>>>--' Testing search formatting  - see test_search.txt '--<<<'
for i in ${out_fmts[@]} ; do echo '>>>--' Testing ${i} '--<<<'
  python ../proggle.py -s "${test_family_ID}" -f ${i} "${test_gff}"
echo "" ; done > test_search.txt
echo '>>>--' Testing search term "(toxin)" '--<<<' >> test_search.txt
python ../proggle.py -s "toxin" -f "Name" "${test_gff}" >> test_search.txt
# 3) search and extract to files
echo -e "\n" '>>>--' Testing search \(writing to files\) - see test_search_tofiles.txt'--<<<'
for i in ${out_fmts[@]} ; do outfile=$(echo "${test_gff}" | sed "s/.gff/_test.${i}/g")
  python ../proggle.py -s "${test_family_ID}" -o "${outfile}" -f ${i} "${test_gff}"
  echo -e "\n" "$(file ${outfile})" "\n" >> test_search_tofiles.txt
  head -30 "${outfile}" >> test_search_tofiles.txt ; rm "${outfile}" ; done
echo -e  '>>>--' Testing feature fasta formatting options '--<<<'
echo "Toxins with 60bp upstream:"
python ../proggle.py -s "toxin" -f "fasta" -us 60 -o "test_toxins_us.fasta" "${test_gff}"
echo "Toxins with 60AA downstream:"
python ../proggle.py -s "toxin" -f "protein" -ds 60 -o "test_toxins_us.aa.fasta" "${test_gff}"
echo "Toxins with 60AA downstream, should overwrite successfully with 100AA line splits:"
python ../proggle.py -s "toxin" -f "protein" -ds 60 -x -o "test_toxins_us.aa.fasta" -fs 100 "${test_gff}"
echo "Repeat of above without overwrite (should raise overwrite error):"
python ../proggle.py -s "toxin" -f "protein" -ds 60 -o "test_toxins_us.aa.fasta" -fs 10 "${test_gff}" # should fail - no overwrite
# 4) testing nearby feature search
# For testing this comprehensively: Test_Genome_00312 is - strand, Test_Genome_00313 is + strand, Test_Genome_00310 is mRNA (no CDS)
echo -e '\n\n>>>--' Testing synteny search features - see test_search.txt '--<<<\n'
echo -e "\n" '>>>--' Testing synteny search features '--<<<' >> test_search.txt
echo -e "\n" '>>>--' Testing search for previous 3 syntenic CDSs '--<<<' >> test_search.txt
python ../proggle.py -n 'before' 3 0 -s "${test_family_ID}" "${test_gff}" >> test_search.txt
echo -e "\n" '>>>--' Testing search for single previous syntenic feature '--<<<' >> test_search.txt
python ../proggle.py -n 'Previous' -s "${test_family_ID}" "${test_gff}" >> test_search.txt
echo -e "\n" '>>>--' Testing search for 3 upstream CDSs '--<<<' >> test_search.txt
python ../proggle.py -n 'upstream' -fh 'ID;product' -s "${test_family_ID}" -us 3 "${test_gff}" >> test_search.txt
echo -e "\n" '>>>--' Testing search for 6th next syntenic CDS, returning ID '--<<<' >> test_search.txt
python ../proggle.py -n 'next' 6 -f ID -s "${test_family_ID}" "${test_gff}" >> test_search.txt
echo -e "\n" '>>>--' Testing search for default no. following features using incorrect but valid input string '--<<<' >> test_search.txt
python ../proggle.py -n 'Accident' -s "${test_family_ID}" "${test_gff}" >> test_search.txt
echo -e "\n" '>>>--' Testing search for neighbouring two flanking syntenic features either side, returning subtab fmt '--<<<' >> test_search.txt
python ../proggle.py -u -n 2 -f 'subtab' -fh 'ID;sequence_length;product' -s "${test_family_ID}" "${test_gff}" >> test_search.txt
echo -e "\n" '>>>--' Testing search for neighbouring two up/downstream syntenic features, ignoring -us and -ds '--<<<' >> test_search.txt
python ../proggle.py -n 'S' 2 2 -us 1000 -ds 2000 -s "${test_family_ID}" "${test_gff}" >> test_search.txt
echo -e "\n" '>>>--' Testing search for 1kbp previous syntenic region and 2kb region following '--<<<' >> test_search.txt
python ../proggle.py -n 'f' 1000 2000 -s "${test_family_ID}" "${test_gff}" >> test_search.txt
echo -e "\n" '>>>--' Testing search for 3kbp up/downstream flanking regions using -us and -ds '--<<<' >> test_search.txt
python ../proggle.py -n 'u' -us 1000 -ds 2000 -s "${test_family_ID}" "${test_gff}" >> test_search.txt
echo -e "\n" '>>>--' Testing search for 1 up/2 downstresam flanking features, returning faa fmt with 9bp us sequence '--<<<' >> test_search.txt
python ../proggle.py -n '~1' 2 -s "${test_family_ID}" -us 3 -f faa -fs 120 "${test_gff}" >> test_search.txt
echo -e '\n\n>>>--' Troubleshooting - checking error in synteny search occurs when using fna output '--<<<\n'
echo -e "\n" '>>>--' Testing search for preceeding flanking region to start of farthest CDS within 3000kb in fna fmt '--<<<' >> test_search.txt
python ../proggle.py -n -3000 -f "fna" -s "${test_family_ID}" "${test_gff}" >> test_search.txt

# Testing editing functionality
# 1) testing GFF editing
echo -e "\n" '>>>--' Testing GFF editing '--<<<'
echo -e "\n" '>>>--' Testing conversion with stat as out_format '--<<<'
python ../proggle.py -f locus_tag "${test_gff}"
echo -e "\n" '>>>--' Testing update and rename contigs - see "${test_gff/.gff/_proggle.gff}" '--<<<'
python ../proggle.py -r -u -o "${test_gff/.gff/_proggle.gff}" -f gff "${test_gff}" 
echo '>>>--' Testing write to directory and without fasta - see test_proggle/ '--<<<'
python ../proggle.py -wf -o "test_proggle/" -f gff "${test_gff}" 
echo '>>>--' Testing write to non-existant subdirectory with contradictory extension - see "test_proggle/subfolder/${test_gff/.gff}.fasta" '--<<<'
python ../proggle.py -o "test_proggle/subfolder/${test_gff/.gff}.fasta" -f faa "${test_gff}" 
# 2) testing adding data files
echo '>>>--' Testing add pangenome data to gff and table - see "${test_gff/.gff/_pan.gff}" and "${test_gff/.gff/_pan.tsv}" '--<<<'
python ../proggle.py -u -a pangenome/gene_presence_absence_panaroo.csv pangenome -o "${test_gff/.gff/_pan.gff}" "${test_gff}"
python ../proggle.py -a pangenome/gene_presence_absence.csv pangenome 1-3,5,8-9 -f table -o "${test_gff/.gff/_pan.tsv}" "${test_gff}"
# 3) file-type conversions
echo -e "\n\n" '>>>--' Testing file type conversions '--<<<'
mkdir test_file_type_conversions
echo "prevent moving gff file itself!" > test_file_type_conversions/"${test_gff}"
for i in ${out_fmts[@]} ; do echo '>>>--' Testing conversion to ${i}'--<<<'
python ../proggle.py -f ${i} "${test_gff}"
mv -nv "${test_gff/.gff}.${i}" test_file_type_conversions/ ; done

# A record of initial testing after 2025/26 NY updates (adding '-c')
# python ../proggle.py -c _3 -f bed -s "aur" "${test_gff}"
# python ../proggle.py -c -3 -f bed -s "aur" "${test_gff}"
# python ../proggle.py -c -1000 -f bed -s "aur" "${test_gff}"
# python ../proggle.py -c 1000 -f bed -s "aur" "${test_gff}"
# python ../proggle.py -c 10 5 -f bed -s "aur" "${test_gff}"
# python ../proggle.py -c '~' 10 5 -f bed -s "aur" "${test_gff}"
# python ../proggle.py -c 's' 10 5 -f bed -s "aur" "${test_gff}"
# python ../proggle.py -c s 10 5 -f bed -s "aur" "${test_gff}"

# Testing contig filtering parameters # now implemented in '-c' as above
# python ../proggle.py -fc 'gnl|UoE|Test_Genome_2' "${test_gff}"
# python ../proggle.py -fc 'gnl|UoE|Test_Genome_2;gnl|UoE|Test_Genome_12' "${test_gff}"
# python ../proggle.py -fc ';gnl|UoE|Test_Genome_2' "${test_gff}"
# python ../proggle.py -fc ';gnl|UoE|Test_Genome_2;gnl|UoE|Test_Genome_12' "${test_gff}"
# python ../proggle.py -fc 10000+ "${test_gff}"
# python ../proggle.py -fc 10000- "${test_gff}"
# python ../proggle.py -fc '21;22;24;26;28' "${test_gff}"
# python ../proggle.py -fc ';21;22;24;26;28' "${test_gff}"

# confirm no changes to raw GFF_feature
diff -qs "${raw_gff}" "${test_gff}"
diff "${raw_gff}" "${test_gff}" > test_gff_changes.txt