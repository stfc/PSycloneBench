
local c = regentlib.c
local cstring = terralib.includec("string.h")

--Terra list containing a list of values we need.
local config_fields_tiles = terralib.newlist({
  --Tilesize options
  {field = "t1", type = int64, default_value = 256, cmd_line = "-t1"},
  {field = "t2", type = int64, default_value = 512, cmd_line = "-t2"}
})

config = terralib.types.newstruct("config")
config.entries:insertall(config_fields_tiles)

-- The following function (get_optional_arg) is under
-- the following copyright and license
-- Copyright 2020 Stanford University
--
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
--
--     http://www.apache.org/licenses/LICENSE-2.0
--
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.
--Taken from pennant_common.rg in the language/examples directory of the Legion runtime.
--This loops through arguments to try to find them, and returns the value if present
local terra get_optional_arg(key : rawstring)
  var args = c.legion_runtime_get_input_args()
  var i = 1
  while i < args.argc do
    if cstring.strcmp(args.argv[i], key) == 0 then
      if i + 1 < args.argc then
        return args.argv[i + 1]
      else
        return nil
      end
    end
    i = i + 1
  end
  return nil
end

terra read_config()

  var conf: config

  --Set defaults - taken from pennant_common.rg
  [config_fields_tiles:map(function(field)
     return quote conf.[field.field] = [field.default_value] end
   end)]

  --For each of the arguments in the config input, get the value using get_optional_arg
  -- and store it in the config structure
  [config_fields_tiles:map(function(field)
      return
      quote 
      var x = get_optional_arg([field.cmd_line])
      if x ~= nil then 
        conf.[field.field] = c.atoll(x) 
      end
      end
  end)]


  return conf
end
