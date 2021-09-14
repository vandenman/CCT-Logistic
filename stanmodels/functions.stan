real get_ordered_logistic(int x, real lt, vector delta) {
  // for debugging purposes
  return ordered_logistic_lpmf(x | lt, delta);
}

int signum(real x) {
  return x > 0 ? 1 :0;
  // if (x > 0) {
  //   return 1;
  // } else {
  //   return -1;
  // }
}
