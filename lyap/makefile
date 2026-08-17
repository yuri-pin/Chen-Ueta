# 1. Compilador e flags de compilacao e debug
FC = gfortran
FFLAGS = -g -fbacktrace -fcheck=all -Wall -fopenmp

#"criar diretório só para escrever menos vezes"
OBJDIR = objeto
BINDIR = executavel

# Flag para o gfortran salvar e buscar os arquivos .mod dentro de "auxiliar"
FFLAGS += -J$(OBJDIR)


# Lista de objetos com o caminho da pasta 'auxiliar'
OBJS_BIF = $(OBJDIR)/func_base_lyap.o $(OBJDIR)/exp_lyap.o


.PHONY: all clean


# 2. Nome do programa executavel final
TARGET = $(BINDIR)/exp_lyap


# 4. Regra padrao do Make
all: $(TARGET)

# Regra genérica para compilar qualquer arquivo .f95 em .o
$(OBJDIR)/%.o: %.f95
	$(FC) $(FFLAGS) -c $< -o $@

# --- PROGRAMA 1: constant_bifurcation ---
$(BINDIR)/exp_lyap: $(OBJS_BIF)
	$(FC) $(FFLAGS) -o  $@ $(OBJS_BIF)


# 7. Dependencia explicita: constant_bifrucation precisa de func_base compilado
$(OBJDIR)/exp_lyap.o: $(OBJDIR)/func_base_lyap.o


# 8. Regra de limpeza
clean:
	rm -rf $(OBJDIR) $(BINDIR)