    def get_MG94Geneconv(self):
        # Change idea here
        # no need to form the huge matrix but only track the non-zero entries
        #Qbasic, distn = self.get_MG94Basic()
        row = []
        col = []
        rate = []
        for i, pair in enumerate(product(self.codon_nonstop, repeat = 2)):
            row_sum = 0
            # use ca, cb, cc to denote codon_a, codon_b, codon_c, where cc != ca, cc != cb
            ca, cb = pair
            sa = self.codon_to_state[ca]
            sb = self.codon_to_state[cb]
            if ca != cb:
                for cc in self.codon_nonstop:
                    if cc == ca or cc == cb:
                        continue
                    sc = self.codon_to_state[cc]
                    # (ca, cb) to (ca, cc)
                    Qbasic = self.get_MG94BasicRate(cb, cc)
                    if Qbasic != 0:
                        row.append((sa, sb))
                        col.append((sa, sc))
                        rate.append(Qbasic)
                        row_sum += Qbasic

                    # (ca, cb) to (cc, cb)
                    Qbasic = self.get_MG94BasicRate(ca, cc)
                    if Qbasic != 0:
                        row.append((sa, sb))
                        col.append((sc, sb))
                        rate.append(Qbasic)
                        row_sum += Qbasic

                        
                # (ca, cb) to (ca, ca)
                row.append((sa, sb))
                col.append((sa, sa))
                Qbasic = self.get_MG94BasicRate(cb, ca)
                if self.isNonsynonymous(cb, ca):
                    Tgeneconv = self.tau * self.omega
                else:
                    Tgeneconv = self.tau
                rate.append(Qbasic + Tgeneconv)
                row_sum += Qbasic + Tgeneconv
                
                # (ca, cb) to (cb, cb)
                row.append((sa, sb))
                col.append((sb, sb))
                Qbasic = self.get_MG94BasicRate(ca, cb)
                rate.append(Qbasic + Tgeneconv)
                row_sum += Qbasic + Tgeneconv

                # Finally add the diagonal
                # (ca, cb) to (ca, cb)
                row.append((sa, sb))
                col.append((sa, sb))
                rate.append(-row_sum)
            else:
                for cc in self.codon_nonstop:
                    if cc == ca:
                        continue
                    sc = self.codon_to_state[cc]

                    # (ca, ca) to (ca,  cc)
                    Qbasic = self.get_MG94BasicRate(ca, cc)
                    if Qbasic != 0:
                        row.append((sa, sb))
                        col.append((sa, sc))
                        rate.append(Qbasic)
                    # (ca, ca) to (cc, ca)
                        row.append((sa, sb))
                        col.append((sc, sa))
                        rate.append(Qbasic)
                        row_sum += 2 * Qbasic
                        
                # Finally add the diagonal
                # (ca, cb) to (ca, cb)
                row.append((sa, sb))
                col.append((sa, sb))
                rate.append(-row_sum)
                        
        process = {'row':row,
                   'col':col,
                   'rate':rate}
        return process

    def get_MG94Sparse(self):
        # Almost the same work flow with the Geneconv one because I am lazy
        row = []
        col = []
        rate = []
        for i, pair in enumerate(product(self.codon_nonstop, repeat = 2)):
            row_sum = 0
            # use ca, cb, cc to denote codon_a, codon_b, codon_c, where cc != ca, cc != cb
            ca, cb = pair
            sa = self.codon_to_state[ca]
            sb = self.codon_to_state[cb]
            if ca != cb:
                for cc in self.codon_nonstop:
                    if cc == ca or cc == cb:
                        continue
                    sc = self.codon_to_state[cc]
                    # (ca, cb) to (ca, cc)
                    Qbasic = self.get_MG94BasicRate(cb, cc)
                    if Qbasic != 0:
                        row.append((sa, sb))
                        col.append((sa, sc))
                        rate.append(Qbasic)
                        row_sum += Qbasic

                    # (ca, cb) to (cc, cb)
                    Qbasic = self.get_MG94BasicRate(ca, cc)
                    if Qbasic != 0:
                        row.append((sa, sb))
                        col.append((sc, sb))
                        rate.append(Qbasic)
                        row_sum += Qbasic

                
                # (ca, cb) to (ca, ca)
                Qbasic = self.get_MG94BasicRate(cb, ca)
                if Qbasic != 0:
                    row.append((sa, sb))
                    col.append((sa, sa))
                    rate.append(Qbasic)
                    row_sum += Qbasic
                
                # (ca, cb) to (cb, cb)
                Qbasic = self.get_MG94BasicRate(ca, cb)
                if Qbasic != 0:
                    row.append((sa, sb))
                    col.append((sb, sb))
                    rate.append(Qbasic)
                    row_sum += Qbasic

                # Finally add the diagonal
                # (ca, cb) to (ca, cb)
                row.append((sa, sb))
                col.append((sa, sb))
                rate.append(-row_sum)
            else:
                for cc in self.codon_nonstop:
                    if cc == ca:
                        continue
                    sc = self.codon_to_state[cc]

                    # (ca, ca) to (ca,  cc)
                    Qbasic = self.get_MG94BasicRate(ca, cc)
                    if Qbasic != 0:
                        row.append((sa, sb))
                        col.append((sa, sc))
                        rate.append(Qbasic)
                    # (ca, ca) to (cc, ca)
                        row.append((sa, sb))
                        col.append((sc, sa))
                        rate.append(Qbasic)

                        row_sum += 2 * Qbasic
                        
                # Finally add the diagonal
                # (ca, cb) to (ca, cb)
                row.append((sa, sb))
                col.append((sa, sb))
                rate.append(-row_sum)
           

        process = {'row':row,
                   'col':col,
                   'rate':rate}
        return process
