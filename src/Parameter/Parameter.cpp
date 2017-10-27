/*
 * Parameter.cpp
 *
 *  Created on: 05/feb/2016
 *      Author: samuele
 */

#include "Parameter.h"

FileParameter::FileParameter() {

}

bool FileParameter::init(string pathFile1, string pathFile2) {
	this->input_files.init(pathFile1, pathFile2);
}

void FileParameter::addSpacedQmer(string nametype, string spaced) {
	this->spaced_seed.push_back(make_pair(nametype, SpacedQmer(spaced)));
}

const PairFiles& FileParameter::getInputFiles() const {
	return input_files;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

VParameter::VParameter() {
	string dir_input = "/media/samuele/Volume/Tesi/input/";
	string dir_sample = dir_input + "Sample/";
	string dir_split = dir_input + "Split/";
	string dir_split_sample = dir_split + "Sample/";
	string dir_metagenoms = "/media/samuele/Volume/Tesi/input/metagenoms/";
	string dir_kraken = dir_metagenoms + "KRAKEN/";
	string dir_kraken_sample = dir_kraken + "Sample/";
	string dir_clark = dir_metagenoms + "CLARK/";
	string dir_clark_sample = dir_clark + "Sample/";
	string dir_real = dir_metagenoms + "SRA/file/";
	string dir_real_sample = dir_real + "Sample/";

//	this->files.push_back(move(FileParameter(dir_sample+"S1.fna.10000")));
//	this->files.push_back(move(FileParameter(dir_split+"S1.fna.1", dir_split+"S1.fna.2")));
//	this->files.push_back(move(FileParameter(dir_input+"R1.fna")));

//	this->files.push_back(move(FileParameter(dir_clark+"A1.10.1000.fq")));
//	this->files.push_back(move(FileParameter(dir_clark+"B1.20.500.fq")));
//	this->files.push_back(move(FileParameter(dir_kraken+"HiSeq_accuracy.fa")));
//	this->files.push_back(move(FileParameter(dir_kraken+"MiSeq_accuracy.fa")));
//	this->files.push_back(move(FileParameter(dir_kraken+"simBA5_accuracy.fa")));
	this->files.push_back(move(FileParameter(dir_kraken+"HiSeq_timing.fa")));
	this->files.push_back(move(FileParameter(dir_kraken+"MiSeq_timing.fa")));
	this->files.push_back(move(FileParameter(dir_kraken+"simBA5_timing.fa")));
//	this->files.push_back(move(FileParameter(dir_kraken_sample+"MiSeq_timing.fa.4000000")));
//	this->files.push_back(move(FileParameter(dir_kraken_sample+"simBA5_timing.fa.3000000")));

//	this->files.push_back(move(FileParameter(dir_split+"S1.fna.1", dir_split+"S1.fna.2")));
//	this->files.push_back(move(FileParameter(dir_split+"S2.fna.1", dir_split+"S2.fna.2")));
//	this->files.push_back(move(FileParameter(dir_split+"S3.fna.1", dir_split+"S3.fna.2")));
//	this->files.push_back(move(FileParameter(dir_split+"S4.fna.1", dir_split+"S4.fna.2")));
//	this->files.push_back(move(FileParameter(dir_split+"S5.fna.1", dir_split+"S5.fna.2")));
	this->files.push_back(move(FileParameter(dir_split+"S6.fna.1", dir_split+"S6.fna.2")));
	this->files.push_back(move(FileParameter(dir_split+"S7.fna.1", dir_split+"S7.fna.2")));
//	this->files.push_back(move(FileParameter(dir_split+"S8.fna.1", dir_split+"S8.fna.2")));
	this->files.push_back(move(FileParameter(dir_split+"S9.fna.1", dir_split+"S9.fna.2")));
	this->files.push_back(move(FileParameter(dir_split+"S10.fna.1", dir_split+"S10.fna.2")));
//	this->files.push_back(move(FileParameter(dir_split_sample+"S10.fna.1500000.1", dir_split_sample+"S10.fna.1500000.2")));

	this->files.push_back(move(FileParameter(dir_split+"MixKraken_abund1.fna.1", dir_split+"MixKraken_abund1.fna.2")));
	this->files.push_back(move(FileParameter(dir_split+"MixKraken_abund2.fna.1", dir_split+"MixKraken_abund2.fna.2")));
//	this->files.push_back(move(FileParameter(dir_split_sample+"MixKraken_abund1.fna.1000000.1", dir_split_sample+"MixKraken_abund1.fna.1000000.2")));
//	this->files.push_back(move(FileParameter(dir_split_sample+"MixKraken_abund2.fna.1000000.1", dir_split_sample+"MixKraken_abund2.fna.1000000.2")));

//	this->files.push_back(move(FileParameter(dir_split+"L1.fna.1", dir_split+"L1.fna.2")));
//	this->files.push_back(move(FileParameter(dir_split+"L2.fna.1", dir_split+"L2.fna.2")));
//	this->files.push_back(move(FileParameter(dir_split+"L3.fna.1", dir_split+"L3.fna.2")));
//	this->files.push_back(move(FileParameter(dir_split+"L4.fna.1", dir_split+"L4.fna.2")));
	this->files.push_back(move(FileParameter(dir_split+"L5.fna.1", dir_split+"L5.fna.2")));
	this->files.push_back(move(FileParameter(dir_split+"L6.fna.1", dir_split+"L6.fna.2")));

//	this->files.push_back(move(FileParameter(dir_input+"R1.fna")));
//	this->files.push_back(move(FileParameter(dir_input+"R2.fna")));
//	this->files.push_back(move(FileParameter(dir_input+"R3.fna")));
//	this->files.push_back(move(FileParameter(dir_input+"R4.fna")));
//	this->files.push_back(move(FileParameter(dir_input+"R5.fna")));
//	this->files.push_back(move(FileParameter(dir_input+"R6.fna")));
	this->files.push_back(move(FileParameter(dir_input+"R7.fna")));
	this->files.push_back(move(FileParameter(dir_input+"R8.fna")));
	this->files.push_back(move(FileParameter(dir_input+"R9.fna")));

	///////////////////////////////////////////////////////
//	this->files.push_back(move(Parameter(directory4+"SRR1804065.fasta.blast_filtered.1", directory4+"SRR1804065.fasta.blast_filtered.2")));
//	this->files.push_back(move(Parameter(directory4+"SRR1804065.fasta.clark_filtered.1", directory4+"SRR1804065.fasta.clark_filtered.2")));
//	this->files.push_back(move(Parameter(directory4+"SRR1804740.fasta.clark_filtered.1", directory4+"SRR1804740.fasta.clark_filtered.2")));
	///////////////////////////////////////////////////////

//	this->spaced_seed.push_back(make_pair("CLARK-S/", SpacedQmer("1111011101110010111001011011111")));
//	this->spaced_seed.push_back(make_pair("CLARK-S/", SpacedQmer("1111101011100101101110011011111")));
//	this->spaced_seed.push_back(make_pair("CLARK-S/", SpacedQmer("1111101001110101101100111011111")));

//	this->spaced_seed.push_back(make_pair("rasbhari_minimizing_overlap_complexity/", SpacedQmer("1111010111010011001110111110111")));
//	this->spaced_seed.push_back(make_pair("rasbhari_minimizing_overlap_complexity/", SpacedQmer("1110111011101111010010110011111")));
//	this->spaced_seed.push_back(make_pair("rasbhari_minimizing_overlap_complexity/", SpacedQmer("1111101001011100111110101101111")));

//	this->spaced_seed.push_back(make_pair("rasbhari_maximizing_sensitivity/", SpacedQmer("1111011110011010111110101011011")));
//	this->spaced_seed.push_back(make_pair("rasbhari_maximizing_sensitivity/", SpacedQmer("1110101011101100110100111111111")));
//	this->spaced_seed.push_back(make_pair("rasbhari_maximizing_sensitivity/", SpacedQmer("1111110101101011100111011001111")));

//	this->spaced_seed.push_back(make_pair("", SpacedQmer("1111111111111111111111")));
//	this->spaced_seed.push_back(make_pair("", SpacedQmer("1010101010101010101010101010101010101010101")));

//	this->spaced_seed.push_back(make_pair("rasbhari_minimizing_overlap_complexity/", SpacedQmer("1111001101110011")));
//	this->spaced_seed.push_back(make_pair("rasbhari_minimizing_overlap_complexity/", SpacedQmer("1100111100011111")));
//	this->spaced_seed.push_back(make_pair("rasbhari_minimizing_overlap_complexity/", SpacedQmer("1111001100101111")));
//	this->spaced_seed.push_back(make_pair("rasbhari_minimizing_overlap_complexity/", SpacedQmer("1110110010110111")));
//	this->spaced_seed.push_back(make_pair("rasbhari_minimizing_overlap_complexity/", SpacedQmer("1101101110111001")));
//	this->spaced_seed.push_back(make_pair("rasbhari_minimizing_overlap_complexity/", SpacedQmer("1011101100101111")));
//	this->spaced_seed.push_back(make_pair("rasbhari_minimizing_overlap_complexity/", SpacedQmer("1111100011010111")));
//	this->spaced_seed.push_back(make_pair("rasbhari_minimizing_overlap_complexity/", SpacedQmer("1101011100011111")));
//	this->spaced_seed.push_back(make_pair("rasbhari_minimizing_overlap_complexity/", SpacedQmer("1111110001011101")));

	this->spaced_seed.push_back(make_pair("rasbhari_minimizing_overlap_complexity/", SpacedQmer("111011100111010001101111001111110011111110111")));
	this->spaced_seed.push_back(make_pair("rasbhari_minimizing_overlap_complexity/", SpacedQmer("101001111001011111011110111111110001010111111")));
	this->spaced_seed.push_back(make_pair("rasbhari_minimizing_overlap_complexity/", SpacedQmer("110101011001100111110101110011001111111111111")));
	this->spaced_seed.push_back(make_pair("rasbhari_minimizing_overlap_complexity/", SpacedQmer("111100011011010010011111111011111111100011111")));
	this->spaced_seed.push_back(make_pair("rasbhari_minimizing_overlap_complexity/", SpacedQmer("100111111111110010010111101111001110110110111")));
	this->spaced_seed.push_back(make_pair("rasbhari_minimizing_overlap_complexity/", SpacedQmer("111101101111100011111110001011101011110111011")));
	this->spaced_seed.push_back(make_pair("rasbhari_minimizing_overlap_complexity/", SpacedQmer("111101001101110011101110101011101110111011111")));
	this->spaced_seed.push_back(make_pair("rasbhari_minimizing_overlap_complexity/", SpacedQmer("110111011111111101101111101010010000011111111")));
	this->spaced_seed.push_back(make_pair("rasbhari_minimizing_overlap_complexity/", SpacedQmer("110100110101111100011111011011111111111110001")));
}
