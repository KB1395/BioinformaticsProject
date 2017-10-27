/*
 * MultiHashFunction.h
 *
 *  Created on: 10/feb/2017
 *      Author: samuele
 */

#ifndef HASH_MULTIHASHFUNCTION_H_
#define HASH_MULTIHASHFUNCTION_H_

#include "HashFunction.h"

inline static void GetHashes_speedup_multi_previous(const string& s_Str,
		const SpacedQmer_Multi& spaced_qmers,
		V_V_Hash_Err& vHashes,
		hash_type (*fConvertion)(char)) {

	auto get_hash = [&](size_t curr_spaced, size_t curr_idx_hash, const PreviusShift_Ext& curr_shift){
		Hash_Err& curr_hash = vHashes.get(curr_spaced, curr_idx_hash);

		if(curr_shift.current_sp_ptr->GetWeight() < curr_shift.GetSize())
			GetHash(s_Str, curr_idx_hash, *curr_shift.current_sp_ptr, curr_hash, fConvertion);
		else
		{
			size_t pos_hash_get = curr_idx_hash-curr_shift.shift_min;//la posizione dell'hash presa è la posizione attuale meno l'indice dello shift dove si fan meno cambiamenti
			const Hash_Err& prev_hash = vHashes.get(curr_shift.prev_qmer, pos_hash_get);

			compute_hash_for_speedup_previous(s_Str,
						curr_shift.current_sp_ptr->GetPosOne(), curr_shift.prev_sp_ptr->GetPosOne(),
						curr_shift,
						prev_hash,
						curr_idx_hash, curr_hash,
						fConvertion);
		}
	};

	//calcolo dimensione v_hash presupponendo che son tutti uguali
	long n_hashes = s_Str.size()-spaced_qmers[0].GetQ()+1;//Presumo lunghezza uguali tra i seed

	vHashes.clear();
	if(n_hashes>0)
	{
		const vector<PreviusShift_Ext_V>& v_shift_min = spaced_qmers.getShiftMin();
		vHashes.resize(spaced_qmers.size(), n_hashes, Hash_Err());

		for(size_t s = 0; s < spaced_qmers.size(); ++s)
		{
			if(!v_shift_min[s][0].isCorrectSpacedPrevious())
				GetHash(s_Str, 0, spaced_qmers[s], vHashes.get(s, 0), fConvertion);//primo da computare a parte
			else
				get_hash(s, 0, v_shift_min[s][0]);
		}

		//TODO: Attenzione a dimensione vettore hash e shift se di dimensioni diverse i seed
		////vHashes[0].size() == vHashes[1].size() == ... == vHashes[n].size()
		////v_shift_min[0].size() == v_shift_min[1].size() == ... == v_shift_min[n].size()
		size_t lim_max = vHashes.size_row(0);
		size_t lim_min = v_shift_min[0].size() < lim_max ? v_shift_min[0].size() : lim_max;

		for(size_t i = 1; i < lim_min; ++i)
			for(size_t ss = 0; ss < spaced_qmers.size(); ++ss)
				get_hash(ss, i, v_shift_min[ss][i]);
		for(size_t i = lim_min; i < lim_max; ++i)//vHashes[0].size() == vHashes[1].size() == ... == vHashes[n].size()
			for(size_t ss = 0; ss < spaced_qmers.size(); ++ss)
				get_hash(ss, i, v_shift_min[ss].back());
	}
	else
		vHashes.resize(spaced_qmers.size(), 0, Hash_Err());
}

inline static void GetHashes_speedup_multi_previous_Rotated(const string& s_Str,
		const SpacedQmer_Multi& spaced_qmers,
		V_V_Hash_Err& vHashes,
		hash_type (*fConvertion)(char)) {

	auto get_hash = [&](size_t curr_spaced, size_t curr_idx_hash, const PreviusShift_Ext& curr_shift){
		Hash_Err& curr_hash = vHashes.get(curr_spaced, curr_idx_hash);

		if(curr_shift.current_sp_ptr->GetWeight() < curr_shift.GetSize())
			GetHash(s_Str, curr_idx_hash, *curr_shift.current_sp_ptr, curr_hash, fConvertion);
		else
		{
			size_t pos_hash_get = curr_idx_hash-curr_shift.shift_min;//la posizione dell'hash presa è la posizione attuale meno l'indice dello shift dove si fan meno cambiamenti
			const Hash_Err& prev_hash = vHashes.get(curr_shift.prev_qmer, pos_hash_get);

			compute_hash_for_speedup_previous(s_Str,
						curr_shift.current_sp_ptr->GetPosOne(), curr_shift.prev_sp_ptr->GetPosOne(),
						curr_shift,
						prev_hash,
						curr_idx_hash, curr_hash,
						fConvertion);
		}
	};

	//calcolo dimensione v_hash presupponendo che son tutti uguali
	long n_hashes = s_Str.size()-spaced_qmers[0].GetQ()+1;//Presumo lunghezza uguali tra i seed

	vHashes.clear();
	if(n_hashes>0)
	{
		const vector<PreviusShift_Ext_V>& v_shift_min = spaced_qmers.getShiftMinRotated();
		vHashes.resize(spaced_qmers.size(), n_hashes, Hash_Err());

		const PreviusShift_Ext_V& zero_prev = v_shift_min[0];
		for(size_t s = 0; s < spaced_qmers.size(); ++s)
		{
			const PreviusShift_Ext& zero = zero_prev[s];
			if(!zero.isCorrectSpacedPrevious())
				GetHash(s_Str, 0, *zero.current_sp_ptr, vHashes.get(s, 0), fConvertion);//primo da computare a parte
			else
				get_hash(s, 0, zero);
		}

		//TODO: Attenzione a dimensione vettore hash e shift se di dimensioni diverse i seed
		////vHashes[0].size() == vHashes[1].size() == ... == vHashes[n].size()
		////v_shift_min[0].size() == v_shift_min[1].size() == ... == v_shift_min[n].size()
		size_t lim_max = vHashes.size_row(0);
		size_t lim_min = v_shift_min.size() < lim_max ? v_shift_min.size() : lim_max;

		for(size_t i = 1; i < lim_min; ++i)
		{
			const PreviusShift_Ext_V& prev_i = v_shift_min[i];
			for(size_t ss = 0; ss < spaced_qmers.size(); ++ss)
				get_hash(ss, i, prev_i[ss]);
		}
		for(size_t i = lim_min; i < lim_max; ++i)//vHashes[0].size() == vHashes[1].size() == ... == vHashes[n].size()
		{
			const PreviusShift_Ext_V& prev_i = v_shift_min.back();
			for(size_t ss = 0; ss < spaced_qmers.size(); ++ss)
				get_hash(ss, i, prev_i[ss]);
		}
	}
	else
		vHashes.resize(spaced_qmers.size(), 0, Hash_Err());
}

#endif /* HASH_MULTIHASHFUNCTION_H_ */
