# g++ -std=c++14 -fpermissive -w -I. で実行 (CommonModule.h/tuple.h はテスト用スタブ)
CXX=g++ -std=c++14 -fpermissive -w -g -I. -I..
all: sweep short
sweep: test_sweep.cpp ../trapezoid.cpp ; $(CXX) $^ -o $@ && ./$@
short: test_short_L.cpp ../trapezoid.cpp ; $(CXX) $^ -o $@ && ./$@
