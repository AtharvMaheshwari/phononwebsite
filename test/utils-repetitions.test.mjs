import assert from 'node:assert/strict';
import { describe, it } from 'node:test';

import { getReasonableRepetitions } from '../src/utils.js';

describe('utils.getReasonableRepetitions', () => {
  it('returns [3,3,3] for tiny systems (<= 10 atoms)', () => {
    assert.deepEqual(getReasonableRepetitions(2), [3, 3, 3]);
    assert.deepEqual(getReasonableRepetitions(10), [3, 3, 3]);
  });

  it('returns [2,2,2] for medium systems (11 to 20 atoms)', () => {
    assert.deepEqual(getReasonableRepetitions(11), [2, 2, 2]);
    assert.deepEqual(getReasonableRepetitions(20), [2, 2, 2]);
  });

  it('returns [1,1,1] for large systems (> 20 atoms)', () => {
    assert.deepEqual(getReasonableRepetitions(21), [1, 1, 1]);
    assert.deepEqual(getReasonableRepetitions(500), [1, 1, 1]);
  });
});
