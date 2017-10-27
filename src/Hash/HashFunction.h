/*
 * HashFunction.h
 *
 *  Created on: 20/lug/2016
 *      Author: samuele
 */

#ifndef HASHFUNCTION_H_
#define HASHFUNCTION_H_

#include "HashType.h"
#include "../Spaced/SpacedQmer_Multi.h"
#include <algorithm>
#include <iostream>

inline static hash_type CharToInt(char ch)
{
	if(ch == 'A')
		return 0;
	if(ch == 'C')
		return 1;
	if(ch == 'G')
		return 2;
	if(ch == 'T')
		return 3;
	return 4; //ERROR CODE
}

inline static hash_type CharToIntComplement(char ch)
{
	if(ch == 'A')
		return 3;
	if(ch == 'C')
		return 2;
	if(ch == 'G')
		return 1;
	if(ch == 'T')
		return 0;
	return 4; //ERROR CODE
}

//Hash per tutti 1 su spaced qmer
inline static void GetHash(const string& s_Str, size_t startQmer, size_t length, Hash_Err& hash_err, hash_type (*fConvertion)(char))
{
	hash_err.reset();
//	#pragma omp parallel for ordered
	for(size_t i = startQmer; i < startQmer + length; ++i)
	{
		hash_type ch = (*fConvertion)(s_Str[i]);
//		#pragma omp ordered
		if(ch == 4) //Errore conversione
			hash_err.push_back_error(i);
		else
			hash_err.hash |= ch << ((i - startQmer) * 2);//OR possibile perchè sommo potenze di 4, OR su posizioni diverse, non c'è riporto
	}
}

//Hash per spaced qmer con *
inline static void GetHash(const string& s_Str, size_t startQmer, const SpacedQmer& spaced_qmer,
		Hash_Err& hash_err, hash_type (*fConvertion)(char))
{
	hash_err.reset();
	const Position& pos_one = spaced_qmer.GetPosOne();
	for(size_t j = 0; j < pos_one.size(); ++j)
	{
		hash_type ch = (*fConvertion)(s_Str[startQmer+pos_one[j]]);
		if(ch == 4) //Errore conversione
			hash_err.push_back_error(j);
		else
			hash_err.hash |= ch << (j * 2);//OR possibile perchè sommo potenze di 4, OR su posizioni diverse, non c'è riporto
	}
}

//Hash veloce con spaced qmer tutti 1
inline static void GetHashes_speedup_previous(const string& s_Str, size_t length,
		Hash_Err_V& vHash, hash_type (*fConvertion)(char)) {
	vHash.clear();
	if(s_Str.size() >= length)
	{
		size_t n_hashes = s_Str.size() - length + 1;
		vHash.resize(n_hashes); //Crea vettore

		GetHash(s_Str, 0, length, vHash[0], fConvertion);//primo da computare a parte
		for(size_t pos=1; pos < vHash.size(); ++pos)
		{
			Hash_Err& prev_hash = vHash[pos-1];
			Hash_Err& curr_hash = vHash[pos];

			//copia hash e sottrai una posizione dal precedente
			curr_hash.hash = prev_hash.hash;
			curr_hash.hash >>= 2; //sposta 2 bit, esce una lettera
			curr_hash.sub_pos_err(1, prev_hash);

			hash_type enter = (*fConvertion)(s_Str[pos+length-1]);
			if(enter == 4)
				curr_hash.push_back_error(length-1);
			else
				curr_hash.hash |= enter << ((length - 1) * 2);	//aggiungi ultimo elemento OR possibile perchè prima ho
																//diviso per 4 e la posizione dove scrivo ha sicuramente 0
		}
	}
}

inline static void GetHashes_naive(const string& s_Str, const SpacedQmer& spaced_qmer,
		Hash_Err_V& vHash, hash_type (*fConvertion)(char))
{
//	bool isAllOne = spaced_qmer.GetWeight() == spaced_qmer.GetQ();
//	if(isAllOne)
//		GetHashes_speedup_previous(s_Str, spaced_qmer.GetQ(), vHash, fConvertion);
//	else
//	{
		vHash.clear();
		if(s_Str.size() >= spaced_qmer.GetQ())
		{
			size_t n_hashes = s_Str.size() - spaced_qmer.GetQ() + 1;
			vHash.resize(n_hashes); //Crea vettore
			#pragma omp parallel for
			for(size_t pos=0; pos < vHash.size(); ++pos)
				GetHash(s_Str, pos, spaced_qmer, vHash[pos], fConvertion);
		}
//	}
}

inline static void compute_hash_for_speedup_previous(const string& s_Str,
													const Position& pos_one_current, const Position& pos_one_prev,
													const PreviousShift& curr_sp_shift,
													const Hash_Err& prev_hash_err,
													size_t idx_curr_hash, Hash_Err& curr_hash_err,
													hash_type (*fConvertion)(char))
{
	//copia hash e errori
	curr_hash_err.hash = prev_hash_err.hash; //Copia hash
	curr_hash_err.hash >>= 2*curr_sp_shift.one_exit;//Shifta correttamente

	//reset one che non fanno più parte dell'hash
	if(!curr_sp_shift.one_to_remove.empty())
	{
		hash_type reset_one = 0;
		for(size_t j = 0; j < curr_sp_shift.one_to_remove.size(); ++j)
			reset_one |= (hash_type)3 << (curr_sp_shift.one_to_remove[j] * 2);
		curr_hash_err.hash &= ~reset_one;
	}

	//Controlla se attualmente hash è corretto
	if(!prev_hash_err.isCorrect())
	{
		long curr_pos_one = 0;
		for(size_t e = 0; e < prev_hash_err.size_error(); ++e)
			if((curr_pos_one = prev_hash_err[e]-curr_sp_shift.one_exit) >= 0)
				if(pos_one_prev[prev_hash_err[e]]-curr_sp_shift.shift_min == pos_one_current[curr_pos_one])
					curr_hash_err.push_back_error(curr_pos_one);//aggiorna posizione errore
	}

	//aggiorna posizioni da cambiare su hash
	for(size_t j = 0; j < curr_sp_shift.one_to_change.size(); ++j)
	{
		const size_t& i_to_change = curr_sp_shift.one_to_change[j];
		size_t index_char = idx_curr_hash+pos_one_current[i_to_change];
		hash_type ch = (*fConvertion)(s_Str[index_char]);
		if(ch == 4) //Errore conversione
			curr_hash_err.push_back_error(i_to_change);
		else
			curr_hash_err.hash |= ch << (i_to_change * 2);//OR possibile perchè sommo potenze di 4, OR su posizioni diverse, non c'è riporto
	}

	//aggiorna rimanenti posizioni da cambiare su hash (quelle uscite son già rimosse)
	//TODO:	si elimina questo pezzo (salta if) se il numero di uno son diversi,
	//in quanto non so dove devo andar ad inserire e rimuovere,
	//NB: l'informazione dove inserire e rimuovere è contenuta tutta
	//sui vettori one_to_change e one_to_remove in quest'ultimo caso
	if(pos_one_current.size() == pos_one_prev.size())
			for(size_t j = pos_one_current.size()-curr_sp_shift.one_exit; j < pos_one_current.size(); ++j)
			{
				size_t index_char = idx_curr_hash+pos_one_current[j];
				hash_type ch = (*fConvertion)(s_Str[index_char]);
				if(ch == 4) //Errore conversione
					curr_hash_err.push_back_error(j);
				else
					curr_hash_err.hash |= ch << (j * 2);//OR possibile perchè sommo potenze di 4, OR su posizioni diverse, non c'è riporto
			}
	////////////////////////////////////////////////////////////////////
	if(!curr_hash_err.isCorrect())
		curr_hash_err.sort_uniq_err();
}

inline static void GetHashes_speedup_previous(const string& s_Str, const SpacedQmer& spaced_qmer,
		Hash_Err_V& vHash, hash_type (*fConvertion)(char)) {
//	bool isAllOne = spaced_qmer.GetWeight() == spaced_qmer.GetQ();
//	if(isAllOne)
//		GetHashes_speedup_previous(s_Str, spaced_qmer.GetQ(), vHash, fConvertion);
//	else
//	{
		auto get_hash = [&](size_t curr_idx_hash, const PreviousShift& curr_shift){
			Hash_Err& curr_hash = vHash[curr_idx_hash];

			if(spaced_qmer.GetWeight() < curr_shift.GetSize())
				GetHash(s_Str, curr_idx_hash, spaced_qmer, curr_hash, fConvertion);
			else
			{
				size_t pos_hash_get = curr_idx_hash-curr_shift.shift_min;//la posizione dell'hash presa è la posizione attuale meno l'indice dello shift dove si fan meno cambiamenti
				const Hash_Err& prev_hash = vHash[pos_hash_get];

				compute_hash_for_speedup_previous(s_Str,
						spaced_qmer.GetPosOne(), spaced_qmer.GetPosOne(),
						curr_shift,
						prev_hash,
						curr_idx_hash, curr_hash,
						fConvertion);
			}
		};

		long n_hashes = s_Str.size()-spaced_qmer.GetQ()+1;
		vHash.clear();
		if(n_hashes>0)
		{
			const V_PreviusShift& shift = spaced_qmer.GetShiftMinChange();
			//Compute hash
			vHash.resize(n_hashes); //Crea vettore

			GetHash(s_Str, 0, spaced_qmer, vHash[0], fConvertion);//primo da computare a parte
			size_t lim_max = vHash.size();
			size_t lim_min = shift.size() < lim_max ? shift.size() : lim_max;
			for(size_t i = 1; i < lim_min; ++i)//Per tutte le posizioni che contemplano gli shift nel primo pezzo di sequenza
				get_hash(i, shift[i]);
			for(size_t i = lim_min; i < lim_max; ++i)
				get_hash(i, shift.back());
		}
//	}
}

#endif /* HASHFUNCTION_H_ */
