//
//  gadgets.hpp
//  gctb
//
//  Created by Jian Zeng on 14/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#ifndef toolbox_hpp
#define toolbox_hpp

#include <sys/time.h>
#include <string>
#include <iostream>

namespace Gadget {

class Timer {
    time_t prev, curr;
public:
    Timer(){
        setTime();
    };
    void setTime(void);
    time_t getTime(void);
    time_t getElapse(void);
    std::string format(const time_t time);
    std::string getDate(void);
    void printElapse(void);
};

}

#endif /* toolbox_hpp */
