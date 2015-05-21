function errors = test_utils()

errors = 0;
errors = errors + test_svdecon();
errors = errors + test_svdsecon();

end

function errors = test_svdecon()

   errors = 0;
   
   x = randn(100,10);
   
   [U,S,V] = svd(x,'econ');   
   [U2,S2,V2] = svdecon(x);
   
   
   if norm(diag(S)-diag(S2),'fro') >1e-9
       errors = errors +1;
      warning('TEST SVDECON 1 pas OK\n');
   else
      fprintf('TEST SVDECON 1 OK\n'); 
   end
   
   [U,S,V] = svd(x','econ');  
   [U2,S2,V2] = svdecon(x');
   
   
  if norm(diag(S)-diag(S2),'fro') >1e-9
       errors = errors +1;
      warning('TEST SVDECON 2 pas OK\n');
   else
      fprintf('TEST SVDECON 2 OK\n'); 
  end
   
    if norm(U2*S2*V2'-x','fro') >1e-9
       errors = errors +1;
      warning('TEST SVDECON 3 pas OK\n');
   else
      fprintf('TEST SVDECON 3 OK\n'); 
   end
end


function errors = test_svdsecon()

   errors = 0;
   
   x = sparse(randn(100,10));
   
   [U,S,V] = svds(x,10);   
   [U2,S2,V2] = svdsecon(x,10);
   
   
   if norm(diag(S)-diag(S2),'fro') >1e-9
       errors = errors +1;
      warning('TEST SVDsECON 1 pas OK\n');
   else
      fprintf('TEST SVDsECON 1 OK\n'); 
   end
   
   [U,S,V] = svds(x',10);  
   [U2,S2,V2] = svdsecon(x',10);
   
   
  if norm(diag(S)-diag(S2),'fro') >1e-9
       errors = errors +1;
      warning('TEST SVDsECON 2 pas OK\n');
   else
      fprintf('TEST SVDsECON 2 OK\n'); 
  end
   
    if norm(U2*S2*V2'-x','fro') >1e-9
       errors = errors +1;
      warning('TEST SVDsECON 3 pas OK\n');
   else
      fprintf('TEST SVDsECON 3 OK\n'); 
   end
end