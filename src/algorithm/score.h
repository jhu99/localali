/**
Author: Jialu Hu
Date: Jun. 25, 2013
File name: algorithm/score.h
Description: Data structure of the scoring function.
**/

#ifndef SCORE_H_
#define SCORE_H_

#include <array>
#include <iostream>

class Score
{
private:
	float _beta;
public:
	std::array<float,5> fscore;
	Score(float);
	~Score(){};
	
};

Score::Score(float beta=1.5)
:_beta(beta)
{
	fscore.fill(0.0);
}

#endif
