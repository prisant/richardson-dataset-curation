-- Does the lab's Top2018 database keep phi/psi of a residue whose neighbour
-- fails the residue filter?  (top2018_mc: mc_filters; top2018_full:
-- allatom_filters.)
--
-- For every residue that passes the filter itself and has a bonded previous
-- and next residue (connectivity.prev_i / next_i), grouped by whether the
-- previous and the next residue pass: how many have phi, psi, and
-- ramachandran_pass_filters = 1.
--
-- Usage: sqlite3 -readonly -header -column top2018_build.db < lab_top2018_neighbour_rule.sql

.print 'mainchain filter (mc_filters)'
SELECT p.mc_pass_filters AS prev_pass, n.mc_pass_filters AS next_pass,
       COUNT(*) AS residues,
       SUM(r.phi IS NOT NULL) AS with_phi,
       SUM(r.psi IS NOT NULL) AS with_psi,
       SUM(r.ramachandran_pass_filters = 1) AS rama_pass
FROM mc_filters f
JOIN connectivity c ON c.chain_i = f.chain_i AND c.residue_i = f.residue_i
JOIN mc_filters p ON p.chain_i = c.chain_i AND p.residue_i = c.prev_i
JOIN mc_filters n ON n.chain_i = c.chain_i AND n.residue_i = c.next_i
JOIN ramachandran r ON r.chain_i = f.chain_i AND r.residue_i = f.residue_i
WHERE f.mc_pass_filters = 1
GROUP BY prev_pass, next_pass;

.print ''
.print 'all-atom filter (allatom_filters, chains in allatom_filters_chainlist)'
SELECT p.all_pass_filters AS prev_pass, n.all_pass_filters AS next_pass,
       COUNT(*) AS residues,
       SUM(r.phi IS NOT NULL) AS with_phi,
       SUM(r.psi IS NOT NULL) AS with_psi,
       SUM(r.ramachandran_pass_filters = 1) AS rama_pass
FROM allatom_filters f
JOIN allatom_filters_chainlist l ON l.chain_i = f.chain_i
JOIN connectivity c ON c.chain_i = f.chain_i AND c.residue_i = f.residue_i
JOIN allatom_filters p ON p.chain_i = c.chain_i AND p.residue_i = c.prev_i
JOIN allatom_filters n ON n.chain_i = c.chain_i AND n.residue_i = c.next_i
JOIN ramachandran r ON r.chain_i = f.chain_i AND r.residue_i = f.residue_i
WHERE f.all_pass_filters = 1
GROUP BY prev_pass, next_pass;

.print ''
.print 'ramachandran_pass_filters vs the mainchain filter of the residue and both neighbours'
.print '(all residues with phi and psi; rule tested: pass = self, prev and next all pass mc)'
SELECT (f.mc_pass_filters = 1 AND p.mc_pass_filters = 1 AND n.mc_pass_filters = 1) AS mc_self_prev_next,
       r.ramachandran_pass_filters AS rama_pass,
       COUNT(*) AS residues
FROM ramachandran r
JOIN mc_filters f ON f.chain_i = r.chain_i AND f.residue_i = r.residue_i
JOIN connectivity c ON c.chain_i = r.chain_i AND c.residue_i = r.residue_i
JOIN mc_filters p ON p.chain_i = c.chain_i AND p.residue_i = c.prev_i
JOIN mc_filters n ON n.chain_i = c.chain_i AND n.residue_i = c.next_i
WHERE r.phi IS NOT NULL AND r.psi IS NOT NULL
GROUP BY mc_self_prev_next, rama_pass;
