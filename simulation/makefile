# windows setup
CXX = g++
CXXFLAGS = -g -Wall
EXEC = simulate
OBJDIR=dist123
OBJS = $(addprefix $(OBJDIR)/, \
	main.o decode.o encode.o PolarCode.o util.o word_error_rate.o build_constraight_matrix.o GaloisFieldPolynomial.o)

${EXEC}: ${OBJS}
	${CXX} ${CXXFLAGS} -o $(OBJDIR)/${EXEC} ${OBJS}


$(OBJDIR)/%.o: %.cpp
	if not exist $(OBJDIR) mkdir $(OBJDIR)
	${CXX} ${CXXFLAGS} -o $@ -c $<

.PHONY: clean

run:
	$(OBJDIR)/${EXEC}

clean:
	del /S /Q $(OBJDIR)
	rmdir $(OBJDIR)