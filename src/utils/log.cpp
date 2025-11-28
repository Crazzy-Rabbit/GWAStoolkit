//
//  log.cpp
//  rsidImpu
//
//  Created by Lulu Shi on 25/11/2025.
//  Copyright © 2025 Lulu Shi. All rights reserved.
//

#include "utils/log.hpp"

// log 文件指针（默认无）
std::ostream* g_log = nullptr;

// 终端输出开关
bool g_log_to_console = true;

// 日志互斥锁
std::mutex g_log_mutex;