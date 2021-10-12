CC:=gcc
src:=rgb_ycc_convert.c
exe:=convert
CFLAGS:=-w -O2

org: $(src)
	$(CC) $(CFLAGS) -DORIGINAL $< -o $(exe)

mp_in: $(src)
	$(CC) $(CFLAGS) -fopenmp -DMP_IN $< -o $(exe)

mp_out: $(src)
	$(CC) $(CFLAGS) -DMP_OUT $< -o $(exe)

.PHONY:clean
clean:
	rm -rf $(exe)

