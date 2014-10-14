require 'Util'

local fieldmt = {}
fieldmt._type = "Field"
fieldmt.__index = fieldmt
local momentummap =  -- momentum lie group for given gauge group
  { GL="GL", SL="TGL", U="AH", SU="TAH", R="R", SR="TR", O="AS", SO="TAS" }

function Field(lattice, opts)
  local f = {}
  f.lattice = lattice
  f.precision = lattice.defaultPrecision
  f.nupdate = 0
  f.wantForce = true
  --f.parents = nil -- hash of [parent]=nupdate
  --f.filename = nil
  --f.metadata = nil
  if opts.kind == "gauge" then
    f.group = lattice.defaultGroup
    f.nc = lattice.defaultNc
  end
  tableCopyTo(f, opts)
  if f.kind == "gauge" then
    if f.momentum then
      f.field = qopqdp.force(f.precision, f.nc, lattice.qdplat)
      f.name = "M" .. f.group .. f.nc
    else
      f.field = qopqdp.gauge(f.precision, f.nc, lattice.qdplat)
      f.name = "GF" .. f.group .. f.nc
    end
  else
    --f.name = "F" .. opts.kind
    error(string.format("unknown Field kind: %s\n", f.kind))
  end
  clearStats(f, "update")
  clearStats(f, "transform")
  return setmetatable(f,fieldmt)
end

function fieldmt:__tostring()
  local s = {}
  tostringRecurse(s, self, {self.lattice})
  if s.done then return "" end
  s[#s+1] = self.name
  s[#s+1] = " = "
  s[#s+1] = self.lattice.name
  s[#s+1] = ":Field {\n"
  s[#s+1] = string.format("  group = \"%s\",\n",self.group)
  s[#s+1] = string.format("  nc = %g\n",self.nc)
  s[#s+1] = "}\n"
  return table.concat(s)
end

function SetTransform(newfs, oldfs, funcs, args)
  for i=1,#newfs do
    local nf = newfs[i]
    if not nf.transform then nf.transform = {} end
    local t = nf.transform
    t.newfs = newfs
    t.oldfs = oldfs
    t.funcs = funcs
    t.args = args
    nf.parents = {}
    for j=1,#oldfs do
      local p = oldfs[j]
      nf.parents[p] = -1
    end
  end
end

-- group, momentum group
function fieldmt:Momentum()
  local fgroup = momentummap[self.group]
  return Field(self.lattice, {kind=self.kind,group=fgroup,nc=self.nc,
			      precision=self.precision,momentum=true})
end

function fieldmt:GetField()
  -- check if parent fields have been updated
  if self.parents then
    local update = false
    for p,n in pairs(self.parents) do
      if p.nupdate > n then
	update = true
	p.nupdate = n
      end
    end
    if update then
      local t = self.transform
      local t0 = clock()
      t.funcs.transform(t.newfs, t.oldfs, t.args)
      updateStats(self, "transform", {seconds=(clock()-t0)})
      self.nupdate = self.nupdate + 1
    end
  end
  return self.field
end

function fieldmt:SetField(f)
  self.field:set(f)
  self.nupdate = self.nupdate + 1
end

function fieldmt:CopyField()
  local f = self:GetField()
  local fc = f:copy()
  return fc
end

function fieldmt:Clone()
  local f = self:GetField()
  local fn = tableCopy(self)
  fn.field = qopqdp.gauge(fn.precision, fn.nc, fn.lattice.qdplat)
  fn.field:set(f)
  return setmetatable(fn,fieldmt)
end

function fieldmt:Set(opts, ...)
  if type(opts)=="table" then
    local t = opts._type
    if t=="Field" then
      self:SetField(opts:GetField())
    elseif t==nil then
      tableCopyTo(self, opts)
    else
      error("unknown type for opts: %s\n", t)
    end
  elseif type(opts)=="string" then -- don't remember what this was for
    self.field[opts](self.field, ...)
    self.nupdate = self.nupdate + 1
  else
    error("unknown type for opts: %s\n", type(opts))
  end
end

function fieldmt:Load(filename)
  local fmd,rmd = self.field:load(filename)
  self.filename = filename
  self.metadata = {file=fmd,record=rmd}
  self.nupdate = self.nupdate + 1
end

function fieldmt:Save(filename, opts)
  if not self.metadata then
    self.metadata = {}
  end
  if not self.metadata.file then
    self.metadata.file =
      "<?xml version=\"1.0\"?>\n<note>generated by qhmc</note>\n"
  end
  if not self.metadata.record then
    self.metadata.record =
      "<?xml version=\"1.0\"?>\n<note>generated by qhmc</note>\n"
  end
  local precision = (opts and opts.precision) or 'F'
  local f = self:GetField()
  f:save(filename,self.metadata.file,self.metadata.record,precision)
  self.filename = filename
end

function fieldmt:Norm2()
  local f = self:GetField()
  return f:norm2()
end

function fieldmt:Random(var) -- ... for variance
  if self.momentum then
    if var then
      self.field:random(var) -- FIXME: only does TAH
    else
      self.field:random() -- FIXME: only does TAH
    end
  else
    self.field:random() -- FIXME: only does SU
  end
  self.nupdate = self.nupdate + 1
end

function fieldmt:Update(momentum, eps)
  local f = self:GetField()
  local m = momentum:GetField()
  local t0 = clock()
  --printf("Field:update: %s\t%g\n", tostring(momentum), eps)
  -- for gauge: exponential, for momentum: add
  f:update(m, eps)
  updateStats(self, "update", {seconds=(clock()-t0)})
  self.nupdate = self.nupdate + 1
end
