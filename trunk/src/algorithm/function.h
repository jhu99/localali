/**
Author: Jialu Hu
Date: Jun. 11, 2013
File name: algorithm/function.h
Description: Searching high-scoring subnetworks.
**/
#pragma once
#ifndef FUNCTION_H_
#define FUNCTION_H_

#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <random>
#include <chrono>
#include <deque>
#include <stack>
#include "verbose.h"
#include "string.h"

/// Convert number types to string type.
template<typename TNum>
std::string
convert_num2str(TNum num)
{
	std::ostringstream buff;
	buff<<num;
	return buff.str();
}

#endif /// FUNCTION_H
