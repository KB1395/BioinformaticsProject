/*
 * FSH.h
 *
 *  Created on: 06/gen/2017
 *      Author: samuele
 */

#ifndef FSH_FSH_H_
#define FSH_FSH_H_

#include "../Hash/MultiHashFunction.h"
#include "../Input/FileScan.h"
#include "../Chrono/Chrono.h"
#include "../Parameter/Parameter.h"

class FSH {
public:
	FSH();
	virtual ~FSH();
	bool load_sequences(const FileParameter& file);
	bool single_run(const SpacedQmer& spaced, bool test_equals);
	void single_save(const FileParameter& file, string dir_output);

	bool multi_run(const vector<SpacedQmer>& multi_spaced, bool test_equals);
	void multi_save(const FileParameter& file, const vector<SpacedQmer>& multi_spaced, string dir_output);

private:
	typedef pair<vector<SpacedQmer>, vector<PreviusShift_Ext_V>> MultiParameter;

	vector<string> to_hash;
	int64_t read_size_max;
	int64_t read_size_min;
	int64_t read_size_mean;

	Chrono::start_end_time times_naive;
	Chrono::start_end_time times_speedup_previous;

	Chrono::start_end_time times_multi_naive;
	Chrono::start_end_time times_multi_speedup_previous;
	Chrono::start_end_time times_multi_speedup_multi_previous;
	bool single_equals = true;

	SpacedQmer_Multi spaced_qmers;
	bool multi_equals = true;

	void single_test_hashes(const SpacedQmer& spaced);
	void single_test_equals(const SpacedQmer& spaced);
	void single_test_naive(const SpacedQmer& spaced);
	void single_test_speedup_previous(const SpacedQmer& spaced);//used previus element

	void multi_test_hashes(const vector<SpacedQmer>& multi_spaced);
	void multi_test_equals(const vector<SpacedQmer>& multi_spaced);
	void multi_test_naive(const vector<SpacedQmer>& multi_spaced);
	void multi_test_speedup_previous(const vector<SpacedQmer>& multi_spaced);//used previus element
	void multi_test_speedup_multi_previous(const vector<SpacedQmer>& multi_spaced);//used previus element of different qmer
};

#endif /* FSH_FSH_H_ */
