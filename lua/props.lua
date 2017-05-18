-- JXY 2016-12-07
-- Combine even-odd and corner to one file
-- Incorporating James's single file writer.

-- ESW 2016-11-14
-- Script file that generates noise sources, dilutes in time and corners, then saves packed traces.
-- Optionally, this script saves the original noise sources and propagators.

require 'common'
require 'smear'
require 'run'

-- Parameters we need to set.
local nx = 48 -- Lattice spatial dimension
local nt = 96 -- Lattice temporal dimension
local mass = 0.0075 -- Light quark mass.
local cfg_num = CONFIG or arg[1] or {printf("NO arg[1]: 0 padded CONFIG_NUMBER\n");os.exit(1)}
local inlat = "./config/config." .. cfg_num ..".lime" -- Input configuration
local outtrace = "./output/trace." .. cfg_num .. ".src" -- Prefix for trace outputs
local outnoise = "./output/noise." .. cfg_num .. ".src" -- Prefix for noise source outputs
local outprop = "./output/prop." .. cfg_num .. ".src" -- Prefix for propagator outputs
local cg_prec = 1e-9 -- Max residual of CG.
local cg_max = 10000 -- Maximum number of CG steps.
local num_stoch = 1; -- How many stochastic sources to generate. (I prefer 6.)
local improved_trace = true; -- Do we use the improved trace estimate or not?
local save_props = true; -- If true, save propagators, if not, only save traces.
local source_type = "Z4"; -- Can also do "Z2", "U1", "Gauss"
local dilute_type = "EO" -- "EO": even odd; "CORNER": corner.  Time dilution always.
local metadata_prefix = "l" .. nx .. ".t" .. nt .. ".m" .. mass .. ".cfg" .. cfg_num; -- Prefix for metadata in saved files.

-----------------
-- BEGIN SETUP --
-----------------


-- These parameters matter!
local prec = 2 -- Use a multi-precision inverter. Default.

-- Start preparing to load the gauge field, spit out basic info.
local latsize = { nx, nx, nx, nt }
local vol = 1
local spatvol = nx*nx*nx;
local seed = seed or os.time()
printf("latsize =")
for k,v in ipairs(latsize) do vol=vol*v; printf(" %i",v); end
printf("\nvolume = %i\n", vol)
printf("mass = %g\n", mass)
printf("seed = %i\n", seed)
printf("cg_prec = %g\n", cg_prec)
printf("num_stoch = %i\n", num_stoch)
printf("source_type = %s\n", source_type)
printf("dilute_type = %s\n", dilute_type)


-- Set up qopqdp.
L = qopqdp.lattice(latsize);
qopqdp.profile(profile or 0);
qopqdp.verbosity(0);
qopqdp.seed(seed);
qopqdp.writeGroupSize(32)   -- default 64

-- Start a timer.
totaltime = qopqdp.dtime()

-- By the way...
local Nc = qopqdp.defaultNc();

function trace_file(i)
  return outtrace .. i .. "." .. source_type .. ".imp" .. (improved_trace and "1" or "0") .. "." .. dilute_type
end
function trace_meta(i)
  return metadata_prefix .. ".type" .. source_type .. ".src" .. i .. ".imp" .. (improved_trace and "1" or "0") .. "." .. dilute_type
end

--------------------
-- SETUP COMPLETE --
--------------------

-------------------------------
-- BEGIN GAUGE FIELD LOADING --
-------------------------------


-- Load the gauge field.
g = qopqdp.gauge();
if inlat then g:load(inlat)
else
  printf("No input configuration specified!\n");
  return 400;
  -- g:random();
end
-- g:unit();

-- Reunitarize, just to be safe.
do
  local devavg,devmax = g:checkSU()
  printf("unitarity deviation avg: %g  max: %g\n", devavg, devmax)
  g:makeSU()
  devavg,devmax = g:checkSU()
  printf("new unitarity deviation avg: %g  max: %g\n", devavg, devmax)
end

-- Print some basic information about the configuration.
function getplaq(g)
  local ps,pt = g:action{plaq=1}
  local lat = qopqdp.lattice()
  local nd,vol = #lat,1
  for i=1,nd do vol=vol*lat[i] end
  local s = 0.25*nd*(nd-1)*vol
  printf("plaq ss: %-.8g  st: %-.8g  tot: %-.8g\n", ps/s, pt/s, 0.5*(ps+pt)/s)
end
getplaq(g);

-- No need to gauge fix!

-----------------------------
-- END GAUGE FIELD LOADING --
-----------------------------

-----------------------------
-- BEGIN ASQTAD SOLVER CFG --
-----------------------------

-- Prepare a smearing setting. This has been verified to be consistent
-- with MILC. Then smear the gauge field!
local smear = {}
smear[#smear+1] = { type="hyp", alpha={0.4,0.5,0.5} }
myprint("smear = ", smear, "\n")

printf("Start smearing.\n");

-- We need to do this because 'smearGauge' expects an action object.
local sg = smearGauge({g = g}, smear);
printf("Smearing done.\n");

-- Set the ASQTAD coefficients. This corresponds to just simple
-- Staggered fermions (supplemented with nHYP smearing)
coeffs = { one_link=1 }

-- Create an asqtad object, set coefficients, etc.
w = qopqdp.asqtad();
w:coeffs(coeffs);
w:set(sg, prec);

-- Set up a function which does an inversion and prints
-- out timing information. Based on actmt.set in asqtadact.lua
function solve_printinfo(w, dest, src, m, res, sub, opts)
  local t0 = qopqdp.dtime();
  w:solve({dest}, src, {m}, res, sub, opts);
  local cgtime = qopqdp.dtime() - t0;
  local flops = w:flops();
  local its = w:its();

  printf("inversion its: %g  secs: %g  Mflops: %g\n", its, cgtime, flops*1e-6/cgtime);
end

---------------------------
-- END ASQTAD SOLVER CFG --
---------------------------

-----------------------------
-- BEGIN GENERATING TRACES --
-----------------------------

printf("Preparing to generate traces.\n");

-- Generate some variables to hold traces, etc, in.

eta = L:colorVector(); -- Hold a noise source.
phi = L:colorVector(); -- Hold a propagator.
trce = L:complex(); -- Hold a trace.
tmpSrc = L:colorVector(); -- Temporary color vector.
tmpSrc2 = L:colorVector(); -- Temporary color vector.
tmpProp = L:colorVector(); -- Another temporary color vector.
tmpTrace = L:complex(); -- Temporary complex.
tmpTrace2 = L:complex(); -- Temporary complex.

for i=1,num_stoch do
  printf("Starting work on noise source %d.\n", i);

  phi:zero();
  trce:zero();

  -- Create the noise source.
  printf("Generating a %s noise source.\n", source_type);
  if (source_type == "Z4") then
    eta:randomZ4();
  elseif (source_type == "Z2") then
    eta:randomZ2();
  elseif (source_type == "U1") then
    eta:randomU1();
  elseif (source_type == "Gauss") then
    eta:random();
    eta:combine({eta}, {1.0/math.sqrt(2.0)}); -- compensate for Gaussian normalization
  else
    printf("ERROR! Invalid noise type %s.\n", source_type);
    io.stdout:flush();

    return 400;
  end

  -- Save the noise source if so desired.
  if (save_props) then
    ftag = i .. "." .. source_type
    pfn = outprop .. ftag
    fmd = metadata_prefix .. ftag
    writer = L:writer(outnoise .. ftag, fmd);
    eta:write(writer, fmd);
    writer:close();
    writer = L:writer(pfn, fmd);
  end

  -- Begin the dilution pattern.
  for t=0,nt-1 do -- Dilute over time.
    tmpSrc:zero();
    tmpSrc:set(eta, "timeslice" .. t);

    local ndilute = 0
    if (dilute_type == "EO") then ndilute = 2
    elseif (dilute_type == "CORNER") then ndilute = 8
    else
      printf("ERROR! Invalid dilute type %s.\n", dilute_type)
      io.stdout:flush();
      return 400;
    end

    for e=0,ndilute-1 do -- Dilute over corners
      local t0 = qopqdp.dtime();
      tmpSrc2:zero();
      local sub = ""
      if (dilute_type == "EO") then sub = e == 0 and "even" or "odd"
      elseif (dilute_type == "CORNER") then sub = "staggered" .. e+8*(t%2)
      else
        printf("ERROR! Invalid dilute type %s.\n", dilute_type)
        io.stdout:flush();
        return 400;
      end

      tmpSrc2:set(tmpSrc, sub);
      printf("Source %d dilution pattern timeslice %d "..dilute_type.." %d.\n", i, t, e);

      phi:zero();
      solve_printinfo(w, phi, tmpSrc2, mass, cg_prec, "all", {prec = prec, restart = cg_max});

      -- Save the propagator if so desired.
      if (save_props) then
        local t0 = qopqdp.dtime();
        --pfn = outprop .. i .. "." .. source_type .. ".timeslice" .. t .. ".corner" .. e
        --fmd = metadata_prefix .. ".type" .. source_type .. ".src" .. i .. .. ".timeslice" .. t .. ".eo" .. e
        --writer = L:writer(pfn, fmd);
        rmd = fmd .. ".timeslice" .. t .. "." .. dilute_type .. e
        phi:write(writer, rmd);
        --writer:close();
        local sptime = qopqdp.dtime() - t0;
        printf("save prop secs: %g\n", sptime);
      end

      -- Perform a contraction.
      tmpTrace:zero();
      if (improved_trace) then
        printf("Computing the improved trace.\n");
        phi:ldot(phi, tmpTrace); -- We rescale by the mass below.
      else
        printf("Computing the unimproved trace.\n");
        tmpSrc2:ldot(phi, tmpTrace);
      end

      -- Update the ldot trace.
      if (improved_trace) then
        trce:combine({trce, tmpTrace}, {1.0, mass});
      else
        trce:combine({trce, tmpTrace}, {1.0, 1.0});
      end

      local ittime = qopqdp.dtime() - t0;
      printf("prop secs: %g\n", ittime);
    end -- for e=0,ndilute-1

  end -- for t=0,nt-1
  if (save_props) then
    writer:close();
  end

  -- Print the per-slice estimate.
  estimates = trce:sum("timeslices");
  for t = 1,nt do
    printf("initsrc %d time %d pbp %.15f\n", i, t, estimates[t].r);
  end

  -- Save the trace.
  writer = L:writer(trace_file(i), trace_meta(i));
  trce:write(writer, trace_meta(i));
  writer:close();

end


---------------------------
-- END GENERATING TRACES --
---------------------------


-------------------------
-- TEST LOADING TRACES --
-------------------------
for i=1,num_stoch do
  printf("Loading trace from source %d.\n", i);

  -- Load the trace.
  reader,file_metadata = L:reader(trace_file(i));
  trace_metadata = trce:read(reader);
  reader:close();

  printf("File metadata for trace %d: %s\n", i, file_metadata);
  printf("Trace metadata for trace %d: %s\n", i, trace_metadata);

  -- Print the per-slice estimate at 0 momentum
  estimates = trce:sum("timeslices");
  for t = 1,nt do
    printf("loadsrc %d mom 0 0 0 time %d pbp %.15f\n", i, t, estimates[t].r);
  end

  -- Print the per-slice estimate at +z momentum
  -- momentum, other field scale, self scale.
  trce:momentum({0,0,1,0},0.0,1.0)
  estimates = trce:sum("timeslices");
  for t = 1,nt do
    printf("loadsrc %d mom 0 0 1 time %d pbp %.15f\n", i, t, estimates[t].r);
  end

end

-----------------------------
-- END TEST LOADING TRACES --
-----------------------------

totaltime = qopqdp.dtime() - totaltime;

printf("Total time: %f seconds.\n", totaltime);

io.stdout:flush()
