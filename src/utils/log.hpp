//
//  log.hpp
//  rsidImpu
//
//  Created by Lulu Shi on 25/11/2025.
//  Copyright © 2025 Lulu Shi. All rights reserved.
//

#ifndef RSIDIMPU_LOG_HPP
#define RSIDIMPU_LOG_HPP

#include <ostream>
#include <iostream>
#include <string>
#include <mutex>

extern std::ostream* g_log;    // 只有在 main.cpp 定义一次
extern bool g_log_to_console;  // 终端输出开关（默认开启）
extern std::mutex g_log_mutex; // 为防止多线程乱序打印，用 mutex

// ------------ logging functions ------------
inline void LOG_INFO(const std::string &msg){
    std::lock_guard<std::mutex> lock(g_log_mutex);

    if (g_log_to_console)
        std::cout << "[INFO] " << msg << std::endl;

    if (g_log) {
        (*g_log) << "[INFO] " << msg << '\n';
    }
}

inline void LOG_WARN(const std::string &msg) {
    std::lock_guard<std::mutex> lock(g_log_mutex);

    if (g_log_to_console)
        std::cout << "[WARN] " << msg << std::endl;

    if (g_log){
        (*g_log) << "[WARN] " << msg << '\n';
    }
}

inline void LOG_ERROR(const std::string &msg) {
    std::lock_guard<std::mutex> lock(g_log_mutex);

    if (g_log_to_console)
        std::cerr << "[ERROR] " << msg << std::endl;
        
    if (g_log) {
        (*g_log) << "[ERROR] " << msg << '\n';
    }
}

#endif