/*
 *  Copyright (C) 2011  Regents of the University of Michigan
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <string>

void testVcfHeader();


const std::string HEADER_LINE="#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003";
const int NUM_SAMPLES = 3;
const std::string SAMPLES[NUM_SAMPLES] = {"NA00001","NA00002","NA00003"};
const int NUM_META_LINES = 18;
const std::string META_LINES[NUM_META_LINES]={
    "##fileformat=VCFv4.1",
    "##fileDate=20090805",
    "##source=myImputationProgramV3.1",
    "##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta",
    "##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=\"Homo sapiens\",taxonomy=x>",
    "##phasing=partial",
    "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">",
    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">",
    "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">",
    "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">",
    "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">",
    "##INFO=<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">",
    "##FILTER=<ID=q10,Description=\"Quality below 10\">",
    "##FILTER=<ID=s50,Description=\"Less than 50% of samples have data\">",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">",
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">",
    "##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">"};

