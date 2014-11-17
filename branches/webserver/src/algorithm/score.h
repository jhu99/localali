/**
Author: Jialu Hu
Date: Jun. 25, 2013
File name: algorithm/score.h
Description: Data structure of the scoring function.
**/
#pragma once

#ifndef SCORE_H_
#define SCORE_H_

#include <array>
#include <iostream>

class Score
{
public:
    /// Evolutionary score for five types of evolutionary events, which include: protein mutation, protein duplication, paralog mutation, interaction deletion, interaction insertion.
	std::array<float,4> fscore;
	Score();
	~Score(){};
	Score& operator+=(Score&);
	Score& operator-=(Score&);
	float sumup();
	void clear();
};

Score::Score()
{
	fscore.fill(0.0);
}

Score& Score::operator+=(Score& another)
{
	for(int i=0;i<4;i++)
	{
		fscore[i]+=another.fscore[i];
	}
	return *this;
}

Score& Score::operator-=(Score& another)
{
	for(int i=0;i<4;i++)
	{
		fscore[i]-=another.fscore[i];
	}
	return *this;
}

float Score::sumup()
{
	float sumscore=0;
	for(int i=0;i<4;i++)
	{
		sumscore+=fscore[i];
	}
	return sumscore;
}

void Score::clear()
{
	fscore.fill(0.0);
}
#endif
