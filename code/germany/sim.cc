void partnerEvent(Agent *a) {
  if (rnd() < (a->partnerSwappiness + a->partner->partnerSwappiness) / 2) { // If a doesn't have a partner, make sure to correct the logic
    find partner for a; // Partner matching algorithm is called here
    if (a->partner is infected) {
      if (rnd() < force_infection) {
        a->infected = true;
      }
    }
  }
}
