# 1. Compilador e flags de compilacao e debug
FC = gfortran
FFLAGS = -g -fbacktrace -fcheck=all -Wall

#"criar diretório só para escrever menos vezes"
OBJDIR = objeto
BINDIR = executavel

# Flag para o gfortran salvar e buscar os arquivos .mod dentro de "auxiliar"
FFLAGS += -J$(OBJDIR)


# Lista de objetos com o caminho da pasta 'auxiliar'
OBJS_BIF = $(OBJDIR)/func_base.o $(OBJDIR)/constant_bifurcation.o
OBJS_IV = $(OBJDIR)/func_base.o $(OBJDIR)/ida_volta.o
OBJS_MLT = $(OBJDIR)/func_base.o $(OBJDIR)/multestabilidade.o
OBJS_RAN = $(OBJDIR)/func_base.o $(OBJDIR)/random_sample.o
OBJS_BASIN = $(OBJDIR)/func_base.o $(OBJDIR)/attraction_basin.o

.PHONY: all clean


# 2. Nome do programa executavel final
TARGET = $(BINDIR)/constant_bifurcation \
         $(BINDIR)/ida_volta \
		 $(BINDIR)/multestabilidade\
		 $(BINDIR)/random_sample\
		 $(BINDIR)/attraction_basin

# 4. Regra padrao do Make
all: $(TARGET)

# Regra genérica para compilar qualquer arquivo .f95 em .o
$(OBJDIR)/%.o: %.f95
	$(FC) $(FFLAGS) -c $< -o $@

# --- PROGRAMA 1: constant_bifurcation ---
$(BINDIR)/constant_bifurcation: $(OBJS_BIF)
	$(FC) $(FFLAGS) -o  $@ $(OBJS_BIF)

# --- PROGRAMA 2: ida_volta ---
$(BINDIR)/ida_volta: $(OBJS_IV)
	$(FC) $(FFLAGS) -o  $@ $(OBJS_IV)

# --- PROGRAMA 3: multestabilidade ---
$(BINDIR)/multestabilidade: $(OBJS_MLT)
	$(FC) $(FFLAGS) -o  $@ $(OBJS_MLT)
	
# --- PROGRAMA 4: random sample ---
$(BINDIR)/random_sample: $(OBJS_RAN)
	$(FC) $(FFLAGS) -o  $@ $(OBJS_RAN)

# --- PROGRAMA 5: attraction basin ---
$(BINDIR)/attraction_basin: $(OBJS_BASIN)
	$(FC) $(FFLAGS) -o  $@ $(OBJS_BASIN)

# 7. Dependencia explicita: constant_bifrucation precisa de func_base compilado
$(OBJDIR)/constant_bifurcation.o: $(OBJDIR)/func_base.o

$(OBJDIR)/ida_volta.o: $(OBJDIR)/func_base.o

$(OBJDIR)/multestabilidade.o: $(OBJDIR)/func_base.o

$(OBJDIR)/random_sample.o: $(OBJDIR)/func_base.o

$(OBJDIR)/attraction_basin.o: $(OBJDIR)/func_base.o

# 8. Regra de limpeza
clean:
	rm -rf $(OBJDIR) $(BINDIR)