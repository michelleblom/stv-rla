# stv-rla
Code for developing RLA configurations for (currently) 2 seat STV elections.

To run with sample size estimation for the generated audits, you will need to
have SHANGRLA (https://github.com/pbstark/SHANGRLA) somewhere in your
PYTHONPATH.

When running any of the scripts, you will see output in the following format:

method used (1Q,2Q,CASEB/C -- represents general method), Instance (path to ballots file), # Candidates, # Voters, Quota, ASN for audit configuration

If the ASN for an instance is 'inf' it means we couldn't find an 'auditable' set of assertions for the instance using the chosen method.


