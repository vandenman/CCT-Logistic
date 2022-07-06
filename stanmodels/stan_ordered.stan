data{
  int<lower = 1> length;
}
parameters{
  ordered[length] x;
}
model {
  x ~ std_normal();
}
