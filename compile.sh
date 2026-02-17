
if [[ $1 == 1 ]]; then
	echo "Compile in debug mode"
	g++ -o penred2Coincidence createCoincidences.cpp common.cpp -O1 -g -fno-omit-frame-pointer -fsanitize=address,bounds,null,leak,undefined
	g++ -o perfect2Coincidence perfectCoin2LM.cpp common.cpp -O1 -g -fno-omit-frame-pointer -fsanitize=address,bounds,null,leak,undefined
	#User export ASAN_OPTIONS="abort_on_error=0:print_stacktrace=1" to run with gdb
else
	echo "Compile in release mode"
	g++ -o penred2Coincidence createCoincidences.cpp common.cpp -O2 -march=native
	g++ -o perfect2Coincidence perfectCoin2LM.cpp common.cpp -O2 -march=native
fi
