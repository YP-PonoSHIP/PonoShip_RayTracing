function out = dbtomag(in)
%This function is used temporarily until constant timeouts are fixed

%   Copyright 2011-2019 The MathWorks, Inc.

%#codegen

coder.extrinsic('db2mag');
out = coder.internal.const(db2mag(in));
