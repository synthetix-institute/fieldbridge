'use strict';
const version=new URL(self.location.href).searchParams.get('physics');
importScripts('model-physics.js'+(version?'?v='+version:''));
self.onmessage=event=>{
  const {id,model,state}=event.data,M=self.FieldBridgeModels;
  try{
    const results=[M.integrate(model,state),M.integrate(model,{...state,initial:state.alternate})];
    if(state.compare)results.push(M.integrate(model,{...state,params:model.params,removed:[]}));
    self.postMessage({id,results});
  }catch(error){self.postMessage({id,error:error.message});}
};
