# AM_RC_CodeReview
Through this study I aim to determine whether different monitoring data streams can be used to better inform Rusty Crayfish dynamics and improve management outcomes.

## The scripts can be divided into a few steps: 
### A. Model Set up 
1. STEP 1: Initial simulation of true abundance= where we simulate true system dynamics (R). We simulate population abundance and data collection (data collection feeds into step 2)
2. STEP 2: First learning process= where we generate our belief of the system (NIMBLE)
3. STEP 3: First decision process= where we select our management and data collection rule. This is based on our belief of the system. 

### B. General Procedure:
4. STEP 4: Simulate Population = simulating true system dynamics
5. STEP 5: Learning process = where we generate our belief of the system using data we collected
6. STEP 6: Update decision process = based on our belief, we update management and data collection rules
