# %%
import os
import io
import matlab.engine

# %% [markdown]
# initializing matlab engine instance

# %%
eng = matlab.engine.start_matlab()
eng.cd("./matlab")
s = eng.genpath("./")

eng.addpath(s)

out = io.StringIO()
err = io.StringIO()

# %% [markdown]
# For some reason running this in parallel uses an ungodly amount of RAM even with a small number of workers, so I'm just opting out of it for now. It doesn't seem to offer the expected increase in throughput either, only around 2x with 16x workers...

# %%
num_seeds = 10240
parallel = 1
correction = 1
try:
    eng.createSeeds(num_seeds, "parallel", parallel, "correction", correction, nargout=0, stdout=out,stderr=err)
finally:
    print(out.getvalue())
    print(err.getvalue())


