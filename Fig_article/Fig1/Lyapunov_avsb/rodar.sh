#!/bin/bash

# Primeiro, compilar todos os arquivos de i=1 a i=100
for i in {1..500}; do
  gcc -O3 lyapunov.c -o i$i.x -lm
  echo
done

# Agora, executar os arquivos em blocos de 7
for i in {1..500}; do
  # Executar os arquivos em blocos de 7 com execução em segundo plano para os primeiros 6, e sequencialmente para o último de cada bloco
  echo $i | ./i$i.x &
  
  # Espera 1 segundo após cada execução
  if (( i % 16 == 0 )); then
    # No final de cada bloco de 7, esperar todas as execuções anteriores terminarem
    wait
  fi
  
  echo
  sleep 1
  rm i$i.x
done
