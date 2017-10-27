/*
 * SpacedQmer.cpp
 *
 *  Created on: 20/lug/2016
 *      Author: samuele
 */

#include "SpacedQmer.h"

SpacedQmer::SpacedQmer()
{
	this->reset("");
}

SpacedQmer::SpacedQmer(string spaced_qmer) {
	this->reset(spaced_qmer);
}

void SpacedQmer::reset(string spaced_qmer) {
	this->spaced_q = spaced_qmer;
	this->SaveIndexOne();
	this->GetShiftMax(this->shift_min_change);
}

void SpacedQmer::SaveIndexOne() {
	this->pos_one.clear();this->pos_one.shrink_to_fit();
	for(size_t i = 0; i < this->spaced_q.length(); ++i)
		if(this->isOne(i))
			this->pos_one.push_back(i);
}

void SpacedQmer::GetShiftMax(V_PreviusShift& shift_max) {
	//reset
	shift_max.clear();
	shift_max.resize(this->spaced_q.size());
	shift_max.shrink_to_fit();

	//calcola posizioni da salvare con i vari shift e poi salva
	size_t init = 0;
	for(size_t i = 1; i < this->spaced_q.size(); ++i)//per tutti gli shift possibili a parte il primo
	{
		//Cerco indice di partenza su vettore uno
		bool find = false;
		for(size_t j = init; j < this->pos_one.size(); ++j)
		{
			if(this->pos_one[j] >= i)
			{
				init = j;
				find = true;
				break;
			}
		}
		if(!find)
			init = this->pos_one.size();//Serve per saltare prossimo ciclo senza altri controlli

		for(size_t j = init; j < this->pos_one.size(); ++j) //per tutte le posizioni del secondo vettore
		{
			//confronta indici incolonnati tenendo conto dello shift, verifica se diversi, se si c'è un operazione da fare in quel punto
			if(this->pos_one[j-init] != this->pos_one[j]-i)
			{
				shift_max[i].one_to_remove.push_back(j-init);
				shift_max[i].one_to_change.push_back(j-init);
			}
		}
//		il prossimo passaggio non è necessario in quanto i bit si spostano tutti di 2*init quindi son già annullati,
//		basta ricordarsi di settarli nella funzione di hash
//		for(size_t j = this->pos_one.size()-init; j < this->pos_one.size(); ++j)//rimanenti da inserire tutti
//			this->shift_min_change[i].one_to_change.push_back(j);

		//Aggiorna minimo e pos_one exit, attenzione con i==1
		shift_max[i].one_exit = init;
		shift_max[i].shift_min = i;//Guess

		const PreviousShift& prev_shift_min = shift_max[shift_max[i-1].shift_min];
		size_t size_previus = prev_shift_min.GetSize();

		size_t size_current = shift_max[i].GetSize();

		if(i > 1 && size_previus < size_current)
			shift_max[i] = shift_max[i-1];
	}
}
