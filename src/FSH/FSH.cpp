/*
 * HashTest.cpp
 *
 *  Created on: 06/gen/2017
 *      Author: samuele
 */

#include "../FSH/FSH.h"

FSH::FSH() {
	// TODO Auto-generated constructor stub

}

FSH::~FSH() {
	// TODO Auto-generated destructor stub
}

bool FSH::load_sequences(const FileParameter& file) {
	this->to_hash.clear();this->to_hash.shrink_to_fit();
	cout << "Load sequences... " << flush;
	bool ret = true;
	FilesScan scans;
	if(scans.init(file.getInputFiles().first.getPath(), file.getInputFiles().second.getPath()))
	{
		auto get_sequence = [](FileScan& filescan, vector<string>& to_hash){
			size_t index_start = to_hash.size();
			to_hash.resize(to_hash.size()+filescan.getSequenceNumber());
			bool ret = true;
			for(size_t i = 0; i < filescan.getSequenceNumber() && ret; ++i)
			{
				Sequence seq;
				if(filescan.getSequenceWithIndex(i, seq))
					to_hash[index_start+i] = seq.getSequence();
				else
					ret = false;
			}
			return ret;
		};
		ret &= get_sequence(scans.first, this->to_hash);
		ret &= get_sequence(scans.second, this->to_hash);
	}
	else
		ret = false;

	//remove not valid
	auto remove_not_valid = [](string& seq){
		for(size_t i = 0; i < seq.size(); ++i)
			if(CharToInt(seq[i]) == 4)
				return true;
		return false;
	};
	auto it = remove_if(this->to_hash.begin(), this->to_hash.end(), remove_not_valid);
	this->to_hash.resize(distance(this->to_hash.begin(), it));

	//count size mean read
	this->read_size_mean = 0;
	for(size_t i = 0; i < this->to_hash.size(); ++i)
		this->read_size_mean = this->read_size_mean + ((int64_t)this->to_hash[i].size() - (int64_t)this->read_size_mean)/(int64_t)(i+1);

	auto min_max = minmax_element(this->to_hash.begin(), this->to_hash.end(), [](string& a, string& b){return a.size()<b.size();});
	this->read_size_min = min_max.first->size();
	this->read_size_max = min_max.second->size();

	cout << "Complete" << endl << flush;
	return ret;
}

bool FSH::single_run(const SpacedQmer& spaced, bool test_equals) {
	this->single_equals = true;
	if(test_equals)
		this->single_test_equals(spaced);
	if(this->single_equals)
		this->single_test_hashes(spaced);
	return this->single_equals;
}

void FSH::single_save(const FileParameter& file, string dir_output) {
	createDirAndSubDir(dir_output);

	auto naive = duration_cast<milliseconds>(this->times_naive.second - this->times_naive.first);
	auto speedup_previous = duration_cast<milliseconds>(this->times_speedup_previous.second - this->times_speedup_previous.first);

	string output_file = dir_output+"times.csv";
	bool exist = file_exist(output_file);
	ofstream output;
	if(exist)
		output.open(output_file, ofstream::out | ofstream::app);
	else
	{
		output.open(output_file, ofstream::out);
		output << "Files,n_seq,m_size,t_naive,t_previous,equals," << endl << flush;
	}
	string str_equals = this->single_equals ? "yes" : "no";
	output << file.getInputFiles().getIdentify() << "," << this->to_hash.size() << "," << this->read_size_mean << "," << naive.count() << "," << speedup_previous.count() << "," << str_equals << "," << endl << flush;
	output.close();
}

void FSH::single_test_hashes(const SpacedQmer& spaced) {
	Chrono chrono;

	cout << "Test naive... " << flush;
	chrono.exe(this, &FSH::single_test_naive, spaced, this->times_naive);
	cout << "Complete" << endl << flush;

	cout << "Test speedup previous... " << flush;
	chrono.exe(this, &FSH::single_test_speedup_previous, spaced, this->times_speedup_previous);
	cout << "Complete" << endl << flush;
}

void FSH::single_test_equals(const SpacedQmer& spaced) {
	cout << "Single test equals hashes... " << flush;
	Hash_Err_V vHash_naive;
	Hash_Err_V vHash_speedup_previous;
	for(size_t i = 0; i < this->to_hash.size() && this->single_equals; ++i)
	{
		GetHashes_naive(this->to_hash[i], spaced, vHash_naive, CharToInt);
		GetHashes_speedup_previous(this->to_hash[i], spaced, vHash_speedup_previous, CharToInt);
		this->single_equals &= (vHash_naive.size() == vHash_speedup_previous.size());
		if(this->single_equals)
		{
			for(size_t j = 0; j < vHash_naive.size(); ++j)
			{
				this->single_equals = (vHash_naive[j].isCorrect() == vHash_speedup_previous[j].isCorrect());
				if(this->single_equals)
				{
					this->single_equals &= (vHash_naive[j].hash == vHash_speedup_previous[j].hash);
					if(!this->single_equals)
						cout<<"";
				}
				else
					cout<<"";
			}
		}
		else
			cout<<"";
	}
	if(this->single_equals)
		cout << "Complete" << endl << flush;
	else
		cout << "Error" << endl << flush;
}

void FSH::single_test_naive(const SpacedQmer& spaced) {
	for(size_t i = 0; i < this->to_hash.size(); ++i)
	{
		Hash_Err_V vHash;
		GetHashes_naive(this->to_hash[i], spaced, vHash, CharToInt);
	}
}

void FSH::single_test_speedup_previous(const SpacedQmer& spaced) {
	for(size_t i = 0; i < this->to_hash.size(); ++i)
	{
		Hash_Err_V vHash;
		GetHashes_speedup_previous(this->to_hash[i], spaced, vHash, CharToInt);
	}
}

bool FSH::multi_run(const vector<SpacedQmer>& multi_spaced, bool test_equals) {
	this->spaced_qmers.init(multi_spaced);

	this->multi_equals = true;
	if(test_equals)
		this->multi_test_equals(multi_spaced);

	if(this->multi_equals)
		this->multi_test_hashes(multi_spaced);
	return this->multi_equals;
}

void FSH::multi_save(const FileParameter& file, const vector<SpacedQmer>& multi_spaced, string dir_output) {
	createDirAndSubDir(dir_output);

	auto naive = duration_cast<milliseconds>(this->times_multi_naive.second - this->times_multi_naive.first);
	auto speedup_previous = duration_cast<milliseconds>(this->times_multi_speedup_previous.second - this->times_multi_speedup_previous.first);
	auto speedup_multi_previous = duration_cast<milliseconds>(this->times_multi_speedup_multi_previous.second - this->times_multi_speedup_multi_previous.first);

	string output_file = dir_output+"multi_times.csv";
	bool exist = file_exist(output_file);
	ofstream output;
	if(exist)
		output.open(output_file, ofstream::out | ofstream::app);
	else
	{
		output.open(output_file, ofstream::out);
		output << "Spaced" << flush;
		for(size_t i = 0; i < multi_spaced.size(); ++i)
			output << "," << multi_spaced[i].toString() << "," <<  endl << flush;

		output << "Files,n_seq,m_size,t_naive,t_previous,t_multi_previous,equals," << endl << flush;
	}
	string str_equals = this->multi_equals ? "yes" : "no";
	output << file.getInputFiles().getIdentify() << "," << this->to_hash.size() << "," << this->read_size_mean << "," << naive.count() << "," << speedup_previous.count() << "," << speedup_multi_previous.count() << "," << str_equals << "," << endl << flush;
	output.close();

	//////changes//////
	output_file = dir_output+"changes.csv";
	output.open(output_file, ofstream::out);
	output << "changes_previous," << endl << flush;
	for(size_t s = 0; s < multi_spaced.size(); ++s)
	{
		const V_PreviusShift& curr_shift = multi_spaced[s].GetShiftMinChange();
		for(size_t i = 0; i < curr_shift.size(); ++i)
			output << curr_shift[i].GetSize() << "," << flush;
		output << endl << flush;
	}
	output << endl << flush;
	output << "changes_multi_previous," << endl << flush;
	const vector<PreviusShift_Ext_V>& v_shift_min = this->spaced_qmers.getShiftMin();
	for(size_t s = 0; s < v_shift_min.size(); ++s)
	{
		for(size_t i = 0; i < v_shift_min[s].size(); ++i)
			output << v_shift_min[s][i].GetSize() << "," << flush;
		output << endl << flush;
	}
	output << endl << flush;
	output << "changes_multi_previous_rotated," << endl << flush;
	const vector<PreviusShift_Ext_V>& v_shift_min_rotated = this->spaced_qmers.getShiftMinRotated();
	for(size_t i = 0; i < v_shift_min_rotated.size(); ++i)
	{
		for(size_t s = 0; s < v_shift_min_rotated[i].size(); ++s)
			output << v_shift_min_rotated[i][s].GetSize() << "," << flush;
		output << endl << flush;
	}
	output.close();
}

void FSH::multi_test_hashes(const vector<SpacedQmer>& multi_spaced) {
	Chrono chrono;

	cout << "Test naive... " << flush;
	chrono.exe(this, &FSH::multi_test_naive, multi_spaced, this->times_multi_naive);
	cout << "Complete" << endl << flush;

	cout << "Test previous... " << flush;
	chrono.exe(this, &FSH::multi_test_speedup_previous, multi_spaced, this->times_multi_speedup_previous);
	cout << "Complete" << endl << flush;

	cout << "Test multi previous... " << flush;
	chrono.exe(this, &FSH::multi_test_speedup_multi_previous, multi_spaced, this->times_multi_speedup_multi_previous);
	cout << "Complete" << endl << flush;
}

void FSH::multi_test_equals(const vector<SpacedQmer>& multi_spaced) {
	cout << "Multi test equals hashes... " << flush;
	vector<vector<Position>> vErr_multi_previous(multi_spaced.size());

	V_V_Hash_Err vHashes_multi_previous;
	Hash_Err_V_V vHashes_single(multi_spaced.size());

	for(size_t i = 0; i < this->to_hash.size() && this->multi_equals; ++i)
	{
		GetHashes_speedup_multi_previous_Rotated(this->to_hash[i], this->spaced_qmers, vHashes_multi_previous, CharToInt);
		for(size_t j = 0; j < multi_spaced.size(); ++j)
			GetHashes_speedup_previous(this->to_hash[i], multi_spaced[j], vHashes_single[j], CharToInt);

		for(size_t j = 0; j < multi_spaced.size(); ++j)
		{
			vector<V_V_Hash_Err::reference_value_type> multi_previous;
			vHashes_multi_previous.get(j, multi_previous);
			this->multi_equals &= multi_previous.size() == vHashes_single[j].size();
			if(!this->multi_equals)
			{
				cout<<"Error size j=" << j <<",prev=" << multi_previous.size() << ",single="<<vHashes_single[j].size()<<endl<<flush;
			}
		}

		if(this->multi_equals)
		{
			for(size_t j = 0; j < multi_spaced.size(); ++j)
			{
				vector<V_V_Hash_Err::reference_value_type> multi_previous;
				vHashes_multi_previous.get(j, multi_previous);
				for(size_t h = 0; h < vHashes_single[j].size() && this->multi_equals; ++h)
				{
					this->multi_equals &= multi_previous[h].get().isCorrect() == vHashes_single[j][h].isCorrect();
					if(this->multi_equals)
					{
						this->multi_equals &= multi_previous[h].get().hash == vHashes_single[j][h].hash;
						if(!this->multi_equals)
						{
							cout<<"Error hash j=" << j <<", h=" << h << ",prev=" <<multi_previous[h].get().hash <<",single="<<vHashes_single[j][h].hash<<endl<<flush;
							cout<<"i="<< i << ","<<this->to_hash[i]<<endl<<flush;
						}
					}
					else
					{
						cout<<"Error charInt j=" << j <<", h=" << h << ",prev=" <<multi_previous[h].get().hash <<",single="<<vHashes_single[j][h].hash<<endl<<flush;
						cout<<"i="<< i << ","<<this->to_hash[i]<<endl<<flush;
					}
				}
			}
		}
		else
			cout<<"";
	}
	if(this->multi_equals)
		cout << "Complete" << endl << flush;
	else
		cout << "Error" << endl << flush;
}

void FSH::multi_test_naive(const vector<SpacedQmer>& multi_spaced) {
	Hash_Err_V_V hash_v(multi_spaced.size());
	for(size_t i = 0; i < this->to_hash.size(); ++i)
		for(size_t j = 0; j < multi_spaced.size(); ++j)
			GetHashes_naive(this->to_hash[i], multi_spaced[j], hash_v[j], CharToInt);
}

void FSH::multi_test_speedup_previous(const vector<SpacedQmer>& multi_spaced) {
	Hash_Err_V_V hash_v(multi_spaced.size());
	for(size_t i = 0; i < this->to_hash.size(); ++i)
		for(size_t j = 0; j < multi_spaced.size(); ++j)
			GetHashes_speedup_previous(this->to_hash[i], multi_spaced[j], hash_v[j], CharToInt);
}

void FSH::multi_test_speedup_multi_previous(const vector<SpacedQmer>& multi_spaced) {
	V_V_Hash_Err hashes_errs;
	hashes_errs.reserve(multi_spaced.size(), this->read_size_max);
	for(size_t i = 0; i < this->to_hash.size(); ++i)
		GetHashes_speedup_multi_previous_Rotated(this->to_hash[i], this->spaced_qmers, hashes_errs, CharToInt);
}

