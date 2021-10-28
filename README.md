# stv-rla
Code for developing RLA configurations for (currently) 2 seat STV elections.

To run you will need to grab assertion_audit_utils.py from Philip Stark's SHANGRLA repository, as we use some functionality from that code, and place it in the same directory as utils.py and audit_general.py.

To get assertion_audit_utils.py, go to the repository: https://github.com/pbstark/SHANGRLA/tree/main/Code

When running any of the scripts, you will see output in the following format:

method used (1Q,2Q,CASEB -- represents general method), Instance (path to ballots file), # Candidates, # Voters, Quota, ASN for audit configuration

If the ASN for an instance is 'inf' it means we couldn't find an 'auditable' set of assertions for the instance using the chosen method.


