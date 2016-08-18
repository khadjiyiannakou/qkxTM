include ../make.inc

all: $(QKXTM)

#%.o: %.cpp $(QKXTM_HDRS)
#	$(CXX) $(INC) $(CXXFLAGS) $< -c -o $@  $(LIBSFLAGS) $(LIBS)   


$(QKXTM): $(QKXTM_OBJS)
	ar cru $@ $(QKXTM_OBJS)

clean:
	-rm -f *.o $(QKXTM)


#%.o: %.c $(QKXTM_HDRS)
#	$(CC) $(INC) $(CFLAGS) $^ -c 