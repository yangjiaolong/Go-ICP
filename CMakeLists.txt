project(GoICP)

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

add_executable(GoICP
	jly_main.cpp
	jly_goicp.cpp
	jly_3ddt.cpp
	matrix.cpp
	ConfigMap.cpp
	StringTokenizer.cpp
	)
