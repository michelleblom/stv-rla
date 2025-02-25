# stv-rla
Code for developing RLA configurations for  STV elections.

To run with sample size estimation for the generated audits, you will need to
have SHANGRLA (https://github.com/pbstark/SHANGRLA) somewhere in your
PYTHONPATH.

Two approaches have been implemented for 2 seat elections:
1) A method (currently called the "1 quota" method) that can be used to generate RLAs for 12 seat STV elections where the first winner wins their seat on first preferences. This means they have more than a quota's worth of votes based on their first preference tally. This is generally quite effective and more often than not we can form a set of assertions that will verify both winners. I have made some improvements to the method compared to what was written in our VOTING'22 paper on the topic, which I am trying to find the time to write up and share.

2) A method (currently called the "general" method) that can be used when the first winner criterion above is not satisfied. This method is going through some flux at the moment. It's generally not effective at forming a set of assertions that will verify both winners. But, I have been looking at what we *can* verify in an audit for these cases. We can often verify the first winner, and beyond this, rule out a number of candidates from being the second winner. The set of assertions being outputted at the moment for these cases needs some refinement, as they can be filtered a bit. I have also realised that the Never Loses (NL) assertions have some additional context that was not clear in the original paper. So, in summary, watch this space!  

An additional variation has been implemented for 2 seat elections:

3) Some US jurisdictions do a batch elimination step before proceeding with the STV algorithm. Basically, this batch eliminates all candidates that cannot mathematically win. We can verify this in an audit. What happens now is that we generate assertions to verify the batch elimination step, and then proceed with either (1) or (2) above.

This repository also contains work for constructing RLAs for 3+ STV elections,
where at least one winner has reached a quota on their first preferences. This
work is in the file audit_n_seats_fwc_simplified.py.	






