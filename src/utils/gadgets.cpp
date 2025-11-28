//
//  gadgets.cpp
//  gctb
//
//  Created by Jian Zeng on 14/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#include "utils/gadgets.hpp"

#include <ctime>
#include <sstream>
#include <iomanip>

void Gadget::Timer::setTime(){
    prev = curr = time(nullptr);
}

time_t Gadget::Timer::getTime(){
    return curr = time(nullptr);
}

time_t Gadget::Timer::getElapse(){
    return curr - prev;
}

std::string Gadget::Timer::format(const time_t t){
    int h = t / 3600;
    int m = (t % 3600) / 60;
    int s = t % 60;

    std::ostringstream oss;
    oss << std::setw(2) << std::setfill('0') << h << ":"
        << std::setw(2) << std::setfill('0') << m << ":"
        << std::setw(2) << std::setfill('0') << s;
    return oss.str();
}

std::string Gadget::Timer::getDate(){
    std::string s = ctime(&curr);
    if (!s.empty() && s.back() == '\n') s.pop_back();
    return s;
}

void Gadget::Timer::printElapse(){
    getTime();
    std::cout << "Time elapse: " << format(getElapse()) << std::endl;
}

